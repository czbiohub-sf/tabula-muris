

library(tidyverse)
library(ggplot2)
library(Seurat)
library(cowplot)
library(dplyr)
# library(viridis,
library(RColorBrewer)
library(here)

droplet_robj = c(
  'droplet_Bladder_seurat_tiss.Robj',
  'droplet_Heart_seurat_tiss.Robj',
  'droplet_Kidney_seurat_tiss.Robj',
  'droplet_Liver_seurat_tiss.Robj',
  'droplet_Lung_seurat_tiss.Robj',
  'droplet_Mammary_seurat_tiss.Robj',
  'droplet_Marrow_seurat_tiss.Robj',
  'droplet_Muscle_seurat_tiss.Robj',
  'droplet_Spleen_seurat_tiss.Robj',
  'droplet_Thymus_seurat_tiss.Robj',
  'droplet_Tongue_seurat_tiss.Robj',
  'droplet_Trachea_seurat_tiss.Robj'
)

facs_robj = c(
  'facs_Aorta_seurat_tiss.Robj',
  'facs_Bladder_seurat_tiss.Robj',
  'facs_Brain_Microglia_seurat_tiss.Robj',
  'facs_Brain_Non-microglia_seurat_tiss.Robj',
  'facs_Colon_seurat_tiss.Robj',
  'facs_Diaphragm_seurat_tiss.Robj',
  'facs_Fat_seurat_tiss.Robj',
  'facs_Heart_seurat_tiss.Robj',
  'facs_Kidney_seurat_tiss.Robj',
  'facs_Liver_seurat_tiss.Robj',
  'facs_Lung_seurat_tiss.Robj',
  'facs_Mammary_seurat_tiss.Robj',
  'facs_Marrow_seurat_tiss.Robj',
  'facs_Muscle_seurat_tiss.Robj',
  'facs_Pancreas_seurat_tiss.Robj',
  'facs_Skin_seurat_tiss.Robj',
  'facs_Spleen_seurat_tiss.Robj',
  'facs_Thymus_seurat_tiss.Robj',
  'facs_Tongue_seurat_tiss.Robj',
  'facs_Trachea_seurat_tiss.Robj'
)

figure_folder = here("20_tissue_tsne_colored_by_expression")


csv = here('00_data_ingest', '16_genes_for_tissue_tsne', 'tSNE_genes_all_tissues.csv')
featureplot_genes = read_csv(csv)
# Unify names
featureplot_genes$Tissue = str_replace_all(featureplot_genes$Tissue, "Brain_microglia", "Brain_Microglia")
featureplot_genes$Tissue = str_replace_all(featureplot_genes$Tissue, "Brain_non_microglia", "Brain_Non-microglia")
unique(featureplot_genes$Tissue)

facs_annotations = read.csv(here('00_data_ingest', '00_facs_raw_data', 'annotations_FACS.csv'), row.names=1)
droplet_annotations = read.csv(here('00_data_ingest', '01_droplet_raw_data', 'annotations_droplets.csv'), row.names=1)
annotations = list(FACS=facs_annotations, Droplet=droplet_annotations)

palette = brewer.pal(3, 'YlGnBu')
cols.use = c(palette[1], palette[3])

feature_plot = function(tiss, tissue_name, platform, genes_to_check) {
  tissue_platform_folder = file.path(figure_folder, tissue_name, platform)
  dir.create(tissue_platform_folder, showWarnings = FALSE, recursive = TRUE)
  print(c('tissue_platform_folder: ', tissue_platform_folder))
  
  for (gene in genes_to_check){
    print(gene)
    # title = sub("_", " ", tissue_name)
    basename = paste0(gene, '.pdf')
    filename = file.path(tissue_platform_folder, basename)
    if( any(gene == tiss@data@Dimnames[[1]]) ){
      print(filename)
      quartz()
      p = FeaturePlot(
        tiss,
        gene,
        pt.size = 1,
        no.legend=FALSE,
        do.return = TRUE,
        no.axes=TRUE,
        cols.use=cols.use,
        )
      # browser()
      # print(paste('p:', p))
      p[[gene]] + guides(color=guide_colorbar(title=""))
      # p + labs(title = title)
      ggsave(filename, width = 2, height = 2)
      dev.off()
    }
  }
}

platforms = list(FACS=facs_robj, Droplet=droplet_robj)
platforms


subset_tissues = c("Heart", "Muscle")

for (platform in names(platforms)){
  print(c("platform: ", platform))
  objects = platforms[[platform]]
  # print(objects)
  
  platform_annotations = annotations[[platform]]
  
  for (robject in objects){
    load(here('00_data_ingest', '10_tissue_robj_downloaded', robject))
    print(c('ls(): ', ls()))
    tissue_of_interest = sub('_seurat_tiss.Robj', '', sub(paste0(tolower(platform), "_"), '', robject))
    print(c("tissue of interest:", tissue_of_interest))
    # if (!startsWith(tissue_of_interest, "Brain")){
    #   next()
    # }
    
    tissue_annotations = platform_annotations[platform_annotations$tissue == tissue_of_interest, ]
    tissue_genes = filter(featureplot_genes, grepl(tissue_of_interest, Tissue, fixed=TRUE))

    genes_to_check = unique(tissue_genes$Gene)
    print(length(genes_to_check))
    if( tissue_of_interest == "Lung" ){
      if( platform == 'FACS'){
        tiss = tissue
      } else {
        tiss = tissue.10X
      }

      # Reassign metadata with cleaned annotations and plot TSNE
      # tiss@meta.data = tissue_annotations
      p = feature_plot(tiss, tissue_of_interest, platform, genes_to_check)
      # rm(list=c('tissue', 'tissue_of_interest', 'tissue.10X'))
    } else{
      # Heart and Aorta annotated some of the same cells but we use the Aorta annotations
      if (any(tissue_of_interest == subset_tissues)){
        tiss = SubsetData(object=tiss, cells.use=rownames(tissue_annotations))
      }

      # Reassign metadata with cleaned annotations and plot TSNE
      # tiss@meta.data = tissue_annotations
      p = feature_plot(tiss, tissue_of_interest, platform, genes_to_check)
      # rm(list=c('tiss', 'tissue_of_interest'))
    }
    
  }
}
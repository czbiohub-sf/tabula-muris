

library(tidyverse)
library(ggplot2)
library(Seurat)
library(cowplot)
library(dplyr)
# library(viridis,
library(RColorBrewer)

plates = c('Aorta/Heart_seurat_tiss.Robj',
            'Bladder/Bladder_seurat_tiss.Robj',
            'Brain_(Non-microglia)/ALL_BRAIN_seurat.Robj',
            "Brain_(Microglia)/microglia_seurat_subtiss.Robj",
            'Colon/Colon_seurat_tiss.Robj',
            "Diaphragm/Muscle_Diaphragm_seurat_tiss.Robj",
            "Fat/Fat_seurat_tiss.Robj",
            'Heart/heart_seurat_tiss.Robj',
            "Kidney/Kidney_seurat_tiss.Robj",
            'Liver/Liver_seurat_tiss-1117.Robj',
            "Lung/SmartSeq2_Lung_seurat_tissue.Robj",
            'Mammary/Mammary_Gland_seurat_tiss.Robj',
            'Marrow/Marrow_seurat_tiss.Robj',
            'Muscle/Muscle_seurat_tiss.Robj',
            'Pancreas/Pancreas_seurat_tiss.Robj',
            'Skin/Skin_seurat_tiss.Robj',
            'Spleen/Spleen_seurat_tiss.Robj',
            'Thymus/Thymus_seurat_tiss.Robj',
            'Tongue/Tongue_seurat_tiss.Robj',
            'Trachea/Trachea_seurat_tiss.Robj'
            )

tenx = c(
  'Bladder/10x_Bladder_seurat_tiss.Robj',
  'Heart/10x_Heart_seurat_tiss.Robj',
  "Kidney/10x_Kidney_seurat_tiss.Robj",
  'Liver/10x_Liver_seurat_3m_201711017.Robj',
  "Lung/10x_Lung_seurat_tissue.10X.Robj",
  'Mammary/10x_Mammary_seurat_tiss.Robj',
  'Marrow/10x_Marrow_seurat_tiss.Robj',
  'Muscle/10x_Muscle_seurat_tiss.Robj',
  'Spleen/10x_Spleen_seurat_tiss.Robj',
  'Thymus/10x_Thymus_final_seurat_tiss.Robj',
  'Tongue/10x_Tongue_seurat_tiss.Robj',
  'Trachea/10x_Trachea_seurat_tiss.Robj'
)

figure_folder = '~/Google Drive/MACA_3mo_manuscript/Main figures/SuppTSNEwithGenes'


featureplot_genes = read_csv('~/code/maca/metadata/tSNE_genes_all_tissues.csv')
# Unify names
featureplot_genes$Tissue = str_replace_all(featureplot_genes$Tissue, "Brain_microglia", "Brain_(Microglia)")
featureplot_genes$Tissue = str_replace_all(featureplot_genes$Tissue, "Brain_non_microglia", "Brain_(Non-microglia)")
unique(featureplot_genes$Tissue)

plates_annotations = read.csv('~/code/maca/metadata/maca_3month_annotations_plates_ontology.csv', row.names=1)
tenx_annotations = read.csv('~/code/maca/metadata/maca_3month_annotations_10x_ontology.csv', row.names=1)
annotations = list(plates=plates_annotations, tenx=tenx_annotations)

palette = brewer.pal(3, 'YlGnBu')
cols.use = c(palette[1], palette[3])

feature_plot = function(tiss, tissue_name, platform, genes_to_check) {
  tissue_platform_folder = file.path(figure_folder, tissue_name, platform)
  dir.create(tissue_platform_folder, showWarnings = FALSE, recursive = TRUE)
  print(c('tissue_platform_folder', tissue_platform_folder))
  
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

platforms = list(FACS=plates, Droplet=tenx)
platforms


setwd("~/Google Drive/MACA uploads/")

object_tissue = c("Lung")
subset_tissues = c("Heart", "Muscle")

for (platform in names(platforms)){
  print(platform)
  objects = platforms[[platform]]
  # print(objects)
  
  platform_annotations = annotations[[platform]]
  
  for (object_name in objects){
    load(object_name)
    # print(ls())
    tissue_of_interest = strsplit(object_name, '/')[[1]][1]
    print(tissue_of_interest)
    if (!startsWith(tissue_of_interest, "Brain")){
      next()
    }
    
    tissue_annotations = platform_annotations[platform_annotations$tissue == tissue_of_interest, ]
    tissue_genes = filter(featureplot_genes, grepl(tissue_of_interest, Tissue, fixed=TRUE))

    genes_to_check = unique(tissue_genes$Gene)
    print(length(genes_to_check))
    if( any(tissue_of_interest == object_tissue)){
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
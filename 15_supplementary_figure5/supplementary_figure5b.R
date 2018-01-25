rm(list=ls())


library(Seurat)
library(cowplot)
library("RColorBrewer")
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


cleaned_annotations = read.csv(here('00_data_ingest', '01_droplet_raw_data', 'annotations_droplets.csv'), row.names=1)


plot_annotated_tsne = function(tiss, object_name, tissue_of_interest) {
  title = sub("_", " ", tissue_of_interest)
  n_annotations = dim(unique(tiss@meta.data['cell_ontology_class']))[1]
  if (n_annotations > 8){
    colors.use = c(brewer.pal(8, 'Set2'), brewer.pal(max(n_annotations-8, 3), 'Dark2'))
  } else{
    colors.use = brewer.pal(max(n_annotations, 3), 'Set2')
  }
  p = TSNEPlot(
    object = tiss,
    do.label = FALSE,
    pt.size = 0.05,
    group.by = 'cell_ontology_class',
    no.legend = TRUE,
    no.axes = TRUE,
    alpha = 0.5,
    do.return = TRUE,
    colors.use=colors.use
  ) #+ geom_point(alpha = 0.1)
  p + labs(title=title)
  ggsave(here('15_supplementary_figure5',
              paste0(
                'figure2b_tsne_',
                tissue_of_interest,
                '.pdf'
              )),
         width = 2,
         height = 2
  )
  return(p)
}

object_tissue = c("Lung")

platform = 'droplet'

for (robject in droplet_robj){
  load(here('00_data_ingest', '10_tissue_robj_downloaded', robject))
  print(ls())
  tissue_of_interest = sub('_seurat_tiss.Robj', '', sub(paste0(platform, "_"), '', robject))
  print(c("tissue of interest:", tissue_of_interest))
  
  tissue_annotations = cleaned_annotations[cleaned_annotations$tissue == tissue_of_interest, ]
  
  if( any(tissue_of_interest == object_tissue)){
    # Reassign metadata with cleaned annotations and plot TSNE
    tissue.10X@meta.data = tissue_annotations
    p = plot_annotated_tsne(tissue.10X, object_name, tissue_of_interest)
    rm(list=c('tissue.10X', 'tissue_of_interest'))
  } else{
    # Reassign metadata with cleaned annotations and plot TSNE
    tiss@meta.data = tissue_annotations
    p = plot_annotated_tsne(tiss, object_name, tissue_of_interest)
    rm(list=c('tiss', 'tissue_of_interest'))
  }
  
}
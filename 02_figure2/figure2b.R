rm(list=ls())



library(Seurat)
library(cowplot)
library("RColorBrewer")
library(here)

facs_robj = c(
  'facs_Aorta_seurat_tiss.Robj',
  'facs_Bladder_seurat_tiss.Robj',
  'facs_Brain_Myeloid_seurat_tiss.Robj',
  'facs_Brain_Non-Myeloid_seurat_tiss.Robj',
  'facs_Large_Intestine_seurat_tiss.Robj',
  'facs_Diaphragm_seurat_tiss.Robj',
  'facs_Fat_seurat_tiss.Robj',
  'facs_Heart_seurat_tiss.Robj',
  'facs_Kidney_seurat_tiss.Robj',
  'facs_Liver_seurat_tiss.Robj',
  'facs_Lung_seurat_tiss.Robj',
  'facs_Mammary_Gland_seurat_tiss.Robj',
  'facs_Marrow_seurat_tiss.Robj',
  'facs_Limb_Muscle_seurat_tiss.Robj',
  'facs_Pancreas_seurat_tiss.Robj',
  'facs_Skin_seurat_tiss.Robj',
  'facs_Spleen_seurat_tiss.Robj',
  'facs_Thymus_seurat_tiss.Robj',
  'facs_Tongue_seurat_tiss.Robj',
  'facs_Trachea_seurat_tiss.Robj'
)

cleaned_annotations = read.csv(here('00_data_ingest', '18_global_annotation_csv', 'annotations_FACS.csv'), row.names=1)


plot_annotated_tsne = function(tiss, object_name, tissue_of_interest) {
  title = sub("_", " ", tissue_of_interest)
  n_annotations = dim(unique(tiss@meta.data['cell_ontology_class']))[1]
  if (n_annotations > 8){
    #colors.use = c(brewer.pal(9, 'Set1'), brewer.pal(max(n_annotations-9, 3), 'Paired'))
    colors.use = c(brewer.pal(8, 'Set2'), brewer.pal(8, 'Dark2'), brewer.pal(8, 'Pastel2'))
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
  ggsave(here('02_figure2',
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

#object_tissue = c("Lung")
#subset_tissues = c("Heart", "Muscle")

platform = 'facs'

for (robject in facs_robj){
  load(here('00_data_ingest', '04_tissue_robj_generated', robject))
  print(ls())
  tissue_of_interest = sub('_seurat_tiss.Robj', '', sub(paste0(platform, "_"), '', robject))
  print(c("tissue of interest:", tissue_of_interest))
  
  tissue_annotations = cleaned_annotations[cleaned_annotations$tissue == tissue_of_interest, ]
  
  # if( tissue_of_interest == "Lung"){
  #   # Reassign metadata with cleaned annotations and plot TSNE
  #   tissue@meta.data = tissue_annotations
  #   p = plot_annotated_tsne(tissue, object_name, tissue_of_interest)
  #   rm(list=c('tissue', 'tissue_of_interest'))
  # } else{
  #   # Heart and Aorta annotated some of the same cells but we use the Aorta annotations
  #   if (any(tissue_of_interest == subset_tissues)){
  #     tiss = SubsetData(object=tiss, cells.use=rownames(tissue_annotations))
  #   }
  # 
  #   # Reassign metadata with cleaned annotations and plot TSNE
  #   tiss@meta.data = tissue_annotations
  #   p = plot_annotated_tsne(tiss, object_name, tissue_of_interest)
  #   rm(list=c('tiss', 'tissue_of_interest'))
  # }
  # 
  p = plot_annotated_tsne(tiss, object_name, tissue_of_interest)
  rm(list=c('tiss', 'tissue_of_interest'))
}

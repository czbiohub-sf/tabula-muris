rm(list = ls())

setwd('~/Google Drive/MACA uploads')

library(Seurat)
library(cowplot)
#
# theme(axis.line=element_blank(),
#       axis.text.x=element_blank(),
#       axis.text.y=element_blank(),
#       axis.ticks=element_blank(),
#       axis.title.x=element_blank(),
#       axis.title.y=element_blank(),
#       legend.position="none",
#       panel.background=element_blank(),
#       panel.border=element_blank(),
#       panel.grid.major=element_blank(),
#       panel.grid.minor=element_blank(),
#       plot.background=element_blank())

objects = c(
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
print(length(objects()))

folder = '~/code/maca/metadata/number_of_cells_reads_genes_10x/'

extract_ngenes_ncells = function(tiss, object) {
  tissue_of_interest = strsplit(object, '/')[[1]][1]
  print(tissue_of_interest)
  
  tissue_metadata = data.frame(c(dim(tiss@scale.data), dim(tiss@raw.data)[2]),
                               row.names = c('n_genes', 'n_cells_pass_qc', 'n_cells_sequenced'))
  colnames(tissue_metadata) = tissue_of_interest
  write.csv(tissue_metadata,
            paste0(folder, tissue_of_interest,
                   '_cell_numbers.csv'))
  
  write.csv(tiss@meta.data[c('nGene', 'nUMI', 'orig.ident')],
            paste0(folder, tissue_of_interest,
                   '_nreads_ngenes.csv'))
}

cleaned_annotations = read.csv('~/code/maca/metadata/maca_3month_annotations_10x_ontology.csv',
                               row.names = 1)

figure_folder = '~/Google Drive/MACA_3mo_manuscript/Main figures/Fig2/10x/'

plot_annotated_tsne = function(tiss, object_name, tissue_of_interest) {
  title = sub("_", " ", tissue_of_interest)
  n_annotations = dim(unique(tiss@meta.data['cell_ontology_class']))[1]
  if (n_annotations > 8){
    colors.use = c(brewer.pal(8, 'Set2'), brewer.pal(n_annotations-8, 'Dark2'))
  } else{
    colors.use = brewer.pal(n_annotations, 'Set2')
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
  ggsave(
    paste0(
      figure_folder,
      'tsne_annotated_',
      tissue_of_interest,
      '.pdf'
    ),
    width = 2,
    height = 2
  )
  return(p)
}

# Lung used a different variable name for their tissue
object_tissue = c("Lung")
#skip_tissues = c("Tongue")

tsne_plots <- list()

i = 0

for (object_name in objects) {
  local({
    # I don't understand why I need to do this but I'm listening to this: https://stackoverflow.com/a/31994539
    i <- i
    tsne_plots <- tsne_plots
  load(object_name)
  print(ls())
  tissue_of_interest = strsplit(object_name, '/')[[1]][1]
  print(tissue_of_interest)
  tissue_annotations = cleaned_annotations[cleaned_annotations$tissue == tissue_of_interest,]
  
  # if (any(tissue_of_interest == skip_tissues)){
  #   next
  # }
  
  if (any(tissue_of_interest == object_tissue)) {
    extract_ngenes_ncells(tissue.10X, object_name)
    
    # Reassign metadata with cleaned annotations and plot TSNE
    tissue.10X@meta.data = tissue_annotations
    p = plot_annotated_tsne(tissue.10X, object_name, tissue_of_interest)
    rm(list = c('tissue.10X', 'tissue_of_interest'))
  } else{
    extract_ngenes_ncells(tiss, object_name)
    
    if (dim(tissue_annotations)[1] == 0) {
      next
    }
    
    # Reassign metadata with cleaned annotations and plot TSNE
    tiss@meta.data = tissue_annotations
    p = plot_annotated_tsne(tiss, object_name, tissue_of_interest)
    rm(list = c('tiss', 'tissue_of_interest'))
  }
  
  # Add to growing list of plots
  i = i + 1
  tsne_plots[[i]] <<- p
  })
  
}
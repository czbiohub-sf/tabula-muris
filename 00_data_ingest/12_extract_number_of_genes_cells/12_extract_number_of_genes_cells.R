library(here)
print(here())

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

platforms = list(facs=facs_robj, droplet=droplet_robj)

extract_ngenes_ncells = function(tiss, object, folder, n_counts='nUMI') {
  suffix = sub(paste0(platform, '_'), '', robject)
  tissue_of_interest = sub('_seurat_tiss.Robj', '', suffix)
  print(tissue_of_interest)
  
  tissue_metadata = data.frame(c(dim(tiss@scale.data), dim(tiss@raw.data)[2]),
                               row.names = c('n_genes', 'n_cells_pass_qc', 'n_cells_sequenced'))
  colnames(tissue_metadata) = tissue_of_interest
  write.csv(tissue_metadata,
            paste0(folder, tissue_of_interest,
                   '_cell_numbers.csv'))
  
  write.csv(tiss@meta.data[c('nGene', n_counts, 'orig.ident')],
            paste0(folder, tissue_of_interest,
                   '_nreads_ngenes.csv'))
}

for (platform in names(platforms)) {
  robjects = platforms[[platform]]
  for (robject in robjects){
    load(here('00_data_ingest', '10_tissue_robj_downloaded', robject))
    print(ls())
    tissue_of_interest = sub('_seurat_tiss.Robj', '', sub(paste0(platform, "_"), '', robject))
    print(c("tissue of interest:", tissue_of_interest))
    
    if (platform == 'facs'){
      folder = paste0(here('00_data_ingest', paste0('13_ngenes_ncells_', platform)), .Platform$file.sep)
    } else {
      folder = paste0(here('00_data_ingest', paste0('14_ngenes_ncells_', platform)), .Platform$file.sep)
    }
    dir.create(folder)
    
    if (platform == 'droplet'){
      n_counts = 'nUMI'
    } else {
      n_counts = 'nReads'
    }
    
    if (tissue_of_interest == 'Lung') {
      # Lung used different naming scheme from rest of tissues
      if (platform == 'droplet'){
        extract_ngenes_ncells(tissue.10X, tissue_of_interest, folder, n_counts=n_counts)
        
        rm(list = c('tissue.10X', 'tissue_of_interest'))
      } else {
        extract_ngenes_ncells(tissue, tissue_of_interest, folder, n_counts=n_counts)
        rm(list=c('tissue', 'tissue_of_interest'))
      }
    } else{
      extract_ngenes_ncells(tiss, tissue_of_interest, folder, n_counts=n_counts)
      rm(list = c('tiss', 'tissue_of_interest'))
    }
  }
}
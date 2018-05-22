library(useful)
library(Seurat)
suppressMessages(library(dplyr))
library(Matrix)
library(ontologyIndex)
suppressMessages(library(tidyverse))

cell_ontology = get_ontology('https://raw.githubusercontent.com/obophenotype/cell-ontology/master/cl-basic.obo', extract_tags='everything')

validate_cell_ontology = function(cell_ontology_class){
  in_cell_ontology = sapply(cell_ontology_class, function(x) is.element(x, cell_ontology$name) || is.na(x))
  if (!all(in_cell_ontology)) {
    message = paste0('"', cell_ontology_class[!in_cell_ontology], '" is not in the cell ontology\n')
    stop(message)
  }
}

convert_to_cell_ontology_id = function(cell_ontology_class){
  return(sapply(cell_ontology_class, function(x) {
    if(is.na(x)){
      x
    }else{
      as.vector(cell_ontology$id[cell_ontology$name == x])[1]
    }
  }))
}

load_tissue_facs = function(tissue_of_interest){
  # Load the per-plate metadata
  plate_metadata_filename = here('00_data_ingest', '00_facs_raw_data', 'metadata_FACS.csv')
  
  plate_metadata <- read.csv(plate_metadata_filename, sep=",", header = TRUE)
  colnames(plate_metadata)[1] <- "plate.barcode"
  
  
  # Load the gene names and set the metadata columns by opening the first file
  filename = here('00_data_ingest', '00_facs_raw_data', 'FACS', paste0(tissue_of_interest, '-counts.csv'))
  
  raw.data = read.csv(filename, sep=",", row.names=1)
  
  plate.barcodes = lapply(colnames(raw.data), function(x) strsplit(strsplit(x, "_")[[1]][1], '.', fixed=TRUE)[[1]][2])
  
  barcode.df = t.data.frame(as.data.frame(plate.barcodes))
  
  rownames(barcode.df) = colnames(raw.data)
  colnames(barcode.df) = c('plate.barcode')
  
  rnames = row.names(barcode.df)
  meta.data <- merge(barcode.df, plate_metadata, by='plate.barcode', sort = F)
  row.names(meta.data) <- rnames
  
  # Sort cells by cell name
  meta.data = meta.data[order(rownames(meta.data)), ]
  raw.data = raw.data[, rownames(meta.data)]
  
  # Find ERCC's, compute the percent ERCC, and drop them from the raw data.
  erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)
  percent.ercc <- Matrix::colSums(raw.data[erccs, ])/Matrix::colSums(raw.data)
  ercc.index <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = FALSE)
  raw.data <- raw.data[-ercc.index,]
  
  # Create the Seurat object with all the data
  tiss <- CreateSeuratObject(raw.data = raw.data, project = tissue_of_interest)
  tiss <- AddMetaData(object = tiss, meta.data)
  tiss <- AddMetaData(object = tiss, percent.ercc, col.name = "percent.ercc")
  
  # Change default name for sums of counts from nUMI to nReads
  colnames(tiss@meta.data)[colnames(tiss@meta.data) == 'nUMI'] <- 'nReads'
  
  # Create metadata columns for cell_ontology_class
  tiss@meta.data[,'free_annotation'] <- NA
  tiss@meta.data[,'cell_ontology_class'] <- NA
  
  #Calculate percent ribosomal genes.
  
  ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = tiss@data), value = TRUE)
  percent.ribo <- Matrix::colSums(tiss@raw.data[ribo.genes, ])/Matrix::colSums(tiss@raw.data)
  tiss <- AddMetaData(object = tiss, metadata = percent.ribo, col.name = "percent.ribo")
  
  tiss <- FilterCells(object = tiss, subset.names = c("nGene", "nReads"), 
                      low.thresholds = c(500, 50000))
  
  tiss <- process_tissue(tiss, 1e6)
  
  return(tiss)
}

stash_annotations = function(tiss, cluster.ids, free_annotation, cell_ontology_class){
  tiss <- StashIdent(object = tiss, save.name = "cluster.ids")

  validate_cell_ontology(cell_ontology_class)
  cell_ontology_id = convert_to_cell_ontology_id(cell_ontology_class)
  
  tiss@meta.data['free_annotation'] <- as.character(plyr::mapvalues(x = tiss@ident, from = cluster.ids, to = free_annotation))
  validate_cell_ontology(cell_ontology_class)
  cell_ontology_id = convert_to_cell_ontology_id(cell_ontology_class)
  
  tiss@meta.data['cell_ontology_class'] <- as.character(plyr::mapvalues(x = tiss@ident, from = cluster.ids, to = cell_ontology_class))
  tiss@meta.data['cell_ontology_id'] <- as.character(plyr::mapvalues(x = tiss@ident, from = cluster.ids, to = cell_ontology_id))
  return(tiss)
}

stash_subtiss_in_tiss = function(tiss, subtiss){
  
  sub.cells = subtiss@cell.names
  
  # Save which cells were in this subset
  subset_cols = sort(grep("subset", colnames(tiss@meta.data), perl=TRUE, value=TRUE))
  if (length(subset_cols) > 0){
    last_subset = subset_cols[length(subset_cols)]
    last_letter = substring(last_subset, 7, 7)
    new_letter = LETTERS[grep(last_letter, LETTERS)+1]
    subset_col = paste0('subset', new_letter)
  } else {
    subset_col = 'subsetA'
  }
  subset_id_col = paste(subset_col, 'cluster.ids', sep='_')

  # Set cells in this subset as TRUE
  tiss@meta.data[, subset_col] = FALSE
  tiss@meta.data[sub.cells, subset_col] = TRUE
  tiss@meta.data[sub.cells, subset_col] = TRUE

  # Save the cluster ids for the subset
  tiss@meta.data[, subset_id_col] = NA
  tiss@meta.data[sub.cells, subset_id_col] = as.numeric(as.vector(subtiss@ident))

  tiss@meta.data[sub.cells, 'free_annotation'] = subtiss@meta.data[,'free_annotation']
  tiss@meta.data[sub.cells, 'cell_ontology_class'] = subtiss@meta.data[,'cell_ontology_class']
  tiss@meta.data[sub.cells, 'cell_ontology_id'] = subtiss@meta.data[,'cell_ontology_id']
  return(tiss)
}


process_tissue = function(tiss, scale){
  tiss <- NormalizeData(object = tiss, scale.factor = scale)
  tiss <- ScaleData(object = tiss)
  tiss <- FindVariableGenes(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5)
  tiss <- RunPCA(object = tiss, do.print = FALSE)
  tiss <- ProjectPCA(object = tiss, do.print = FALSE)
}

load_tissue_droplet = function(tissue_of_interest){
    
  # read the metadata to get the channels we want
  droplet_metadata_filename = here('00_data_ingest', '01_droplet_raw_data', 'metadata_droplet.csv')
  
  droplet_metadata <- read.csv(droplet_metadata_filename, sep=",", header = TRUE)
  colnames(droplet_metadata)[1] <- "channel"

  tissue_metadata = filter(droplet_metadata, tissue == tissue_of_interest)[,c('channel','tissue','subtissue','mouse.sex', 'mouse.id')]
  
  subfolder = paste0(tissue_of_interest, '-', tissue_metadata$channel[1])
  raw.data <- Read10X(data.dir = here('00_data_ingest', '01_droplet_raw_data', 'droplet', subfolder))
  colnames(raw.data) <- lapply(colnames(raw.data), function(x) paste0(tissue_metadata$channel[1], '_', x))
  meta.data = data.frame(row.names = colnames(raw.data))
  meta.data['channel'] = tissue_metadata$channel[1]
  
  if (length(tissue_metadata$channel) > 1){
    # Some tissues, like Thymus and Heart had only one channel
    for(i in 2:nrow(tissue_metadata)){
      subfolder = paste0(tissue_of_interest, '-', tissue_metadata$channel[i])
      new.data <- Read10X(data.dir = here('00_data_ingest', '01_droplet_raw_data', 'droplet', subfolder))
      colnames(new.data) <- lapply(colnames(new.data), function(x) paste0(tissue_metadata$channel[i], '_', x))
      
      new.metadata = data.frame(row.names = colnames(new.data))
      new.metadata['channel'] = tissue_metadata$channel[i]
      
      raw.data = cbind(raw.data, new.data)
      meta.data = rbind(meta.data, new.metadata)
    }
  }
  
  rnames = row.names(meta.data)
  meta.data <- merge(meta.data, tissue_metadata, sort = F)
  row.names(meta.data) <- rnames

  # Order the cells alphabetically to ensure consistency.

  ordered_cell_names = order(colnames(raw.data))
  raw.data = raw.data[,ordered_cell_names]
  meta.data = meta.data[ordered_cell_names,]

  # Find ERCC's, compute the percent ERCC, and drop them from the raw data.
  erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)
  percent.ercc <- Matrix::colSums(raw.data[erccs, ])/Matrix::colSums(raw.data)
  ercc.index <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = FALSE)
  raw.data <- raw.data[-ercc.index,]
  
  # Create the Seurat object with all the data
  tiss <- CreateSeuratObject(raw.data = raw.data, project = tissue_of_interest)
  
  tiss <- AddMetaData(object = tiss, meta.data)
  tiss <- AddMetaData(object = tiss, percent.ercc, col.name = "percent.ercc")

  ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = tiss@data), value = TRUE)
  percent.ribo <- Matrix::colSums(tiss@raw.data[ribo.genes, ])/Matrix::colSums(tiss@raw.data)
  tiss <- AddMetaData(object = tiss, metadata = percent.ribo, col.name = "percent.ribo")

  # Create metadata columns for annotations
  tiss@meta.data[,'free_annotation'] <- NA
  tiss@meta.data[,'cell_ontology_class'] <- NA
  
  tiss <- FilterCells(object = tiss, subset.names = c("nGene", "nUMI"), 
                      low.thresholds = c(500, 1000))
  tiss <- process_tissue(tiss, 1e4)
  
  return(tiss)
}


save_annotation_csv = function(tiss, tissue_of_interest, method='facs', additional_cols=c()){
  if (method == "facs"){
    batch_name_column = 'plate.barcode'
  } else {
    batch_name_column = 'channel'
  }
  tiss@meta.data['cell'] = rownames(tiss@meta.data)
  
  filename = here('00_data_ingest', '03_tissue_annotation_csv',
                    paste0(tissue_of_interest, "_", method, "_annotation.csv"))

  # Set columns commonly used across all methods
  common_cols = c('cell', 'tissue', 'subtissue', 'cell_ontology_class',
                'cell_ontology_id', 'free_annotation', 'cluster.ids', 'mouse.sex',
                'mouse.id', 'tSNE_1', 'tSNE_2')


  # Get any columns that say "subset"
  subset_cols = sort(grep("subset", colnames(tiss@meta.data), perl=TRUE, value=TRUE))

  columns = c(common_cols, batch_name_column, subset_cols, additional_cols)

  write_csv(FetchData(tiss, columns), filename)
            
}

compare_previous_annotation = function(tiss, tissue_of_interest, method='facs', filename=NULL){
  if (is.null(filename)){
    filename = here('00_data_ingest', '03_tissue_annotation_csv', 
                    paste0(tissue_of_interest, "_", method, "_annotation.csv"))
  }
  if (file.exists(filename)){
    previous_annotation = read.csv(filename, stringsAsFactors = FALSE)
    cols = c('free_annotation', 'cell_ontology_class')
    for (col in cols){
      previous_col = paste0('previous_', col)
      tiss@meta.data[, previous_col] <- "NA"
      tiss@meta.data[as.character(previous_annotation$X), previous_col] <- previous_annotation[, col]
      print(table(tiss@meta.data[, previous_col]))
      print(table(tiss@meta.data[, previous_col], tiss@ident))

    }
  }
  return(tiss)
}

previous_annotation_table = function(tiss, method='facs'){
  previous_cols = c('previous_free_annotation', 'previous_cell_ontology_class')
  for (previous_col in previous_cols){
    print(previous_col)
    print(table(tiss@meta.data[, previous_col], tiss@ident))
  }
}

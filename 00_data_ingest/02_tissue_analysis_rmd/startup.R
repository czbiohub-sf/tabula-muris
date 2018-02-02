library(here)
library(useful)
library(Seurat)
library(dplyr)
library(Matrix)
library(ontologyIndex)
cell_ontology = get_ontology('https://raw.githubusercontent.com/obophenotype/cell-ontology/master/cl-basic.obo', extract_tags='everything')

validate_cell_ontology = function(cell_ontology_class){
  in_cell_ontology = sapply(cell_ontology_class, function(x) is.element(x, cell_ontology$name))
  if (!all(in_cell_ontology)) {
    message = paste0('"', cell_ontology_class[!in_cell_ontology], '" is not in the cell ontology\n')
    stop(message)
  }
}
convert_to_cell_ontology_id = function(cell_ontology_class){
  return(sapply(cell_ontology_class, function(x) as.vector(cell_ontology$id[cell_ontology$name == x])[1]))
}

save_dir = here('00_data_ingest', 'tissue_robj')
# 22_markers

This folder hosts the code and csv files used to compare the differentially expressed genes per cell type across Tabula muris and microwell-seq.

## Available files

### Code
* `FindMarkers.Rmd`: this function re-runs `FindAllMarkers` but instead of comparing across clusters compares across cell types. 
* `methods_comparison.m`: this function finds the overlapping differentially expressed genes across cell types
* `vennX.m`: auxiliary function, draws proportional Venn diagrams. Adapted from Jeremy Heil's [vennX](https://www.mathworks.com/matlabcentral/fileexchange/6116-proportional-venn-diagrams) 

### Cell type gene signature
#### Tabula muris
* `method_organ_cell_ontology_class_classes.csv`
* `method_organ_cell_ontology_class_markers.csv`

#### Mouse cell atlas (microwell-seq)
* `organ_Han.csv`

### Methods comparison results
* `organ_cell_ontologies.csv`: cell types that we matched across the three technologies
* `organ_genes.csv`: gene lists corresponding to all sections of the Venn diagrams per cell type

## Getting started

### From tissue "robj"

After you get each `robj` per organ per method in folder `04_tissue_robj_generated` you can run `FindMarkers.Rmd` to re-run Seurat's `FindAllMarkers` function and compute the differentially expressed genes between the different cell types annottated for both methods used in Tabula muris.

This outputs the csv files with the **Cell type gene signatures** for each tissue and method.

### From cell type gene signature csv files

If using the available csv files, you can directly run the function `methods_comparison.m`. It will ask the user to choose which tissue is of interest and will output the corresponding Venn diagram showing the overlap among the differentially expressed genes per cell type across the three methodologies. This information is also written into the **Methods comparison results** csv files.

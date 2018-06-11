## ----setup---------------------------------------------------------------
library(knitr)
knit_hooks$set(optipng = hook_optipng)

## ------------------------------------------------------------------------
library(here)
source(here('30_tissue_supplement_figures', 'supplemental_figures.R'))
save_folder = here('30_tissue_supplement_figures', 'Liver', 'facs')
dir.create(save_folder, recursive=TRUE)
method = "facs"

tissue_of_interest = 'Liver'
filename = paste0('facs_',tissue_of_interest, '_seurat_tiss.Robj')
load(here('00_data_ingest', '04_tissue_robj_generated', filename))

# Make sure cluster ids are numeric
tiss@meta.data[, 'cluster.ids'] = as.numeric(tiss@meta.data[, 'cluster.ids'])

# Concatenate original cell ontology class to free annotation
cell_ontology_class = tiss@meta.data$cell_ontology_class
cell_ontology_class[is.na(cell_ontology_class)] = "NA"

free_annotation = sapply(tiss@meta.data$free_annotation,
    function(x) { if (is.na(x)) {return('')} else return(paste(":", x))},
    USE.NAMES = FALSE)
tiss@meta.data[, "free_annotation"] = paste(cell_ontology_class,
    free_annotation, sep='')

additional.group.bys = sort(c())

group.bys = c(standard.group.bys, additional.group.bys)

genes_to_check = c("Alb", "Apoa1", "Ass1", "Cd19", "Cd68", "Cd74", "Cd79a", "Cd79b", "Cdh1", "Clec4f", "Cxcr6", "Cyp2e1", "Cyp2f2", "Emr1", "Glul", "Gstp1", "Gulo", "Gzma", "Hal", "Hamp", "Il2rb", "Irf7", "Kdr", "Nkg7", "Nrp1", "Oat", "Oit3", "Pck1", "Pecam1", "Serpina1c", "Ttr", "Ubb", "Zap70")

## ----use-optipng, optipng='-o7'------------------------------------------
dot_tsne_ridge(tiss, genes_to_check, save_folder, prefix = prefix,
    group.bys = group.bys, method = method)

## ------------------------------------------------------------------------
#tiss.markers <- FindAllMarkers(object = tiss, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
#filename = file.path(save_folder, paste(prefix, 'findallmarkers.csv', sep='_'))
#write.csv(tiss.markers, filename)


## ----setup---------------------------------------------------------------
library(knitr)
knit_hooks$set(optipng = hook_optipng)

## ------------------------------------------------------------------------
library(here)
source(here('30_tissue_supplement_figures', 'supplemental_figures.R'))
save_folder = here('30_tissue_supplement_figures', 'Mammary_Gland', 'facs')
dir.create(save_folder, recursive=TRUE)
method = "facs"

tissue_of_interest = 'Mammary_Gland'
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

genes_to_check = c("Aldh1a3", "Ccl5", "Cd14", "Cd19", "Cd3e", "Cd3g", "Cd55", "Cd74", "Cd79a", "Cdh5", "Csf1r", "Elf5", "Esam", "Esr1", "Fn1", "Krt14", "Krt17", "Krt18", "Krt19", "Krt5", "Krt8", "Lgr5", "Ly6c1", "Nkg7", "Pecam1", "Pgr", "Prlr", "Ptprc", "Vim")

## ----use-optipng, optipng='-o7'------------------------------------------
dot_tsne_ridge(tiss, genes_to_check, save_folder, prefix = prefix,
    group.bys = group.bys, method = method)

## ------------------------------------------------------------------------
#tiss.markers <- FindAllMarkers(object = tiss, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
#filename = file.path(save_folder, paste(prefix, 'findallmarkers.csv', sep='_'))
#write.csv(tiss.markers, filename)

## ----optipng='-o7'-------------------------------------------------------
in_SubsetA = tiss@meta.data$cell_ontology_class == "basal cell"
in_SubsetA[is.na(in_SubsetA)] = FALSE


## ----optipng='-o7'-------------------------------------------------------
SubsetA.cells.use = tiss@cell.names[in_SubsetA]
write(paste("Number of cells in SubsetA subset:", length(SubsetA.cells.use)), stderr())
SubsetA.n.pcs = 10
SubsetA.res.use = 0.5
SubsetA.perplexity = 30
SubsetA.genes_to_check = c("Areg", "Cited1", "Csn3", "Fos", "Id1", "Id2", "Igfbp2", "Igfpb4", "Prlr", "Procr", "Vcam1")
SubsetA.group.bys = c(group.bys, "subsetA_cluster.ids")
SubsetA.tiss = SubsetData(tiss, cells.use=SubsetA.cells.use, )
SubsetA.tiss <- SubsetA.tiss %>% ScaleData() %>% 
  FindVariableGenes(do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5) %>%
  RunPCA(do.print = FALSE)
SubsetA.tiss <- SubsetA.tiss %>% FindClusters(reduction.type = "pca", dims.use = 1:SubsetA.n.pcs, 
    resolution = SubsetA.res.use, print.output = 0, save.SNN = TRUE) %>%
    RunTSNE(dims.use = 1:SubsetA.n.pcs, seed.use = 10, perplexity=SubsetA.perplexity)


## ----optipng='-o7'-------------------------------------------------------
colors.use = c('LightGray', 'Coral')
tiss@meta.data[, "SubsetA"] = "(Not in subset)"
tiss@meta.data[SubsetA.tiss@cell.names, "SubsetA"] = "SubsetA" 
filename = make_filename(save_folder, prefix="SubsetA", 'highlighted', 
    'tsneplot_allcells')
p = TSNEPlot(
  object = tiss,
  do.return = TRUE,
  group.by = "SubsetA",
  no.axes = TRUE,
  pt.size = 1,
  no.legend = TRUE,
  colors.use = colors.use
) + coord_fixed(ratio = 1) +
    xlab("tSNE 1") + ylab("tSNE 2")
ggsave(filename, width = 4, height = 4)

filename = make_filename(save_folder, prefix="SubsetA", 'highlighted', 
    'tsneplot_allcells_legend')
# Plot TSNE again just to steal the legend
p = TSNEPlot(
    object = tiss,
    do.return = TRUE,
    group.by = "SubsetA",
    no.axes = TRUE,
    pt.size = 1,
    no.legend = FALSE,
    label.size = 8,
    colors.use = colors.use
    ) + coord_fixed(ratio = 1) +
    xlab("tSNE 1") + ylab("tSNE 2")

# Initialize an empty canvas!
ggdraw()
# Draw only the legend
ggdraw(g_legend(p))
ggsave(filename, width = 8, height = 4)
dev.off()


## ----optipng='-o7'-------------------------------------------------------
dot_tsne_ridge(SubsetA.tiss, SubsetA.genes_to_check,
    save_folder, prefix = "SubsetA-Basal", group.bys = SubsetA.group.bys, 
    "facs")



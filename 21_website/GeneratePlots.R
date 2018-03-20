## ------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(here)
library(parallel)
library(foreach)
library(pryr)

## ------------------------------------------------------------------------
source(here("00_data_ingest", "02_tissue_analysis_rmd", "boilerplate.R"))

## ------------------------------------------------------------------------
genes = read_csv("gene_list.csv", col_names = c("gene_names"))$gene_names
metas = c("cell_ontology_class", "free_annotation", "mouse.sex")
meta_display_names = c("Cell Ontology Class", "Free Annotation", "Sex")

number_of_cores = 16

## ------------------------------------------------------------------------
plotgene_tsne <- function(tiss, gene, tissue, method){
  gene_R_safe = paste0("`",gene,"`")
  if(method == "facs"){
    legend_name = "ln(1+CPM)"
  }  else{
    legend_name = "ln(1+CP10k)"
  }
  lims = FetchData(tiss, c('tSNE_1', 'tSNE_2', gene)) %>% summarize(xmin = min(tSNE_1), xmax = max(tSNE_1), ymin = min(tSNE_2), ymax = max(tSNE_2))
  plot_min = min(lims$xmin, lims$ymin)
  plot_max = min(lims$xmax, lims$ymax)
  FetchData(tiss, c('tSNE_1', 'tSNE_2', gene)) %>% ggplot(aes_string(x = 'tSNE_1', y = 'tSNE_2', color = gene_R_safe)) +
    geom_point(size = 0.5) +
    scale_colour_gradient(low = "lightgrey", high = "blue", name = legend_name) +
    xlim(plot_min, plot_max) + ylim(plot_min, plot_max) + coord_fixed(ratio = 1) +
    xlab("tSNE 1") + ylab("tSNE 2")

  ggsave(here('21_website','images', paste0(tissue, "-", method, "-", gene,"-tsne", ".png")))
}

plotgene_vln <- function(tiss, gene, tissue, method){
  if(method == "facs"){
    legend_name = "ln(1+CPM)"
  }  else{
    legend_name = "ln(1+CP10k)"
  }
  VlnPlot(tiss, gene, group.by = 'cell_ontology_class') +
  xlab("Cell Ontology Class") + ylab(paste0("Expression: ", legend_name)) + ggtitle("") +
    coord_flip()
  ggsave(here('21_website','images', paste0(tissue, "-", method, "-", gene,"-vln", ".png")), width = 14, height = 7)
}

plotmeta <- function(tiss, tissue, method, index){
  meta = metas[index]
  legend_name = meta_display_names[index]

  lims = FetchData(tiss, c('tSNE_1', 'tSNE_2', gene)) %>% summarize(xmin = min(tSNE_1), xmax = max(tSNE_1), ymin = min(tSNE_2), ymax = max(tSNE_2))
  plot_min = min(lims$xmin, lims$ymin)
  plot_max = min(lims$xmax, lims$ymax)

  meta_R_safe = paste0("`",meta,"`")
  FetchData(tiss, c('tSNE_1', 'tSNE_2', meta)) %>% ggplot(aes_string(x = 'tSNE_1', y = 'tSNE_2', color = meta_R_safe)) +
    geom_point(size = 0.5) +
    xlim(plot_min, plot_max) + ylim(plot_min, plot_max) + coord_fixed(ratio = 1) +
    scale_colour_discrete(name = legend_name) +
    xlab("tSNE 1") + ylab("tSNE 2")

  ggsave(here('21_website',"images", paste0(tissue, "-", method, "-", meta, "-tsne.png")))
}

## ------------------------------------------------------------------------
tissue_plots <- function(tissue, method){
  load(here("00_data_ingest","04_tissue_robj_generated", paste0(method, "_", tissue, "_", "seurat_tiss.Robj")))

  gene_only_plotgene_tsne = pryr::partial(plotgene_tsne, tiss=tiss, tissue=tissue, method=method)
  gene_only_plotgene_vln = pryr::partial(plotgene_vln, tiss=tiss, tissue=tissue, method=method)
  meta_only_plotmeta = pryr::partial(plotmeta, tiss=tiss, tissue=tissue, method=method)

  mclapply(genes, gene_only_plotgene_tsne, mc.preschedule=TRUE, mc.cores=number_of_cores)
  mclapply(genes, gene_only_plotgene_vln, mc.preschedule=TRUE, mc.cores=number_of_cores)
  mclapply(1:length(metas), meta_only_plotmeta, mc.preschedule=TRUE, mc.cores=number_of_cores)
}

## ------------------------------------------------------------------------
ptm <- proc.time()

tissue_plots("Pancreas", "facs")
tissue_plots("Muscle", "droplet")
tissue_plots("Mammary", "facs")

proc.time() - ptm


# Title     : Generate TSNE, dotplots, and violinplots for tissue supplement
# Objective : make lots of plots
# Created by: olgabot
# Created on: 3/19/18

# Load packages
library(useful)
library(Seurat)
library(dplyr)
library(Matrix)
library(ontologyIndex)
library(tidyverse)
library(cowplot)
library(grid)

# Commonly used variables by all Rmd files
standard.group.bys = c("cell_ontology_class", "free_annotation", 'cluster.ids')
prefix = 'allcells'

CHUNKSIZE = 12

# Extract legend from a ggplot
# Stolen from https://gist.github.com/crsh/be88be19233f1df4542aca900501f0fb
# Who stole it from: http://stackoverflow.com/a/12041779/914024
gglegend <- function(x){
  tmp <- ggplot_gtable(ggplot_build(x))
  leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
  tmp$grobs[[leg]]
}

# Make all supplemental figures for a tiss object
dot_tsne_violin = function(tiss, genes_to_check, save_folder, prefix, group.bys){
    for (group.by in group.bys ){
      write(paste('group.by:', group.by, '   prefix:', prefix), stderr())

      # If this column is all NAs, then skip it
      if (all(is.na(tiss@meta.data[, group.by]))){
        next
      }

      # TSNE only
      filename = file.path(save_folder, paste(prefix, group.by,
        'tsneplot.pdf', sep='_'))
      p = TSNEPlot(object = tiss, do.return = TRUE, group.by = group.by,
        no.axes=TRUE, pt.size=1, no.legend=TRUE)
      ggsave(filename, width = 2, height = 2)
      dev.off()

      # Legend for TSNE
      filename = file.path(save_folder, paste(prefix, group.by,
        'tsneplot_legend.pdf', sep='_'))
      p = TSNEPlot(object = tiss, do.return = TRUE, group.by = group.by,
        no.axes=TRUE, pt.size=1, no.legend=FALSE)
      # dev.off()
      grid.draw(gglegend(p))
      ggsave(filename, width = 2, height = 2)
      dev.off()

      chunked_genes = split(genes_to_check, ceiling(seq_along(genes_to_check)/CHUNKSIZE))
      n_chunks = length(chunked_genes)
      for (name in names(chunked_genes)){
        write(paste('\tchunk:', name), stderr())

        genes = chunked_genes[[name]]
        # Dotplot - enrichment of gene expression in group.by with dot size
        filename = file.path(save_folder, paste(prefix, group.by,
          paste0('dotplot_', name, '-of-', n_chunks, '.pdf'), sep='_'))
      p = DotPlot(tiss, genes, col.max = 2.5, plot.legend = T,
        do.return = T, group.by = group.by) #+ coord_flip()
      ggsave(filename, width = 11, height = 8)
      dev.off()

      # Violinplot - enrichment of gene expression in group.by with double wide histogram
      filename = file.path(save_folder, paste(prefix, group.by,
          paste0('violinplot_', name, '-of-', n_chunks, '.pdf'), sep='_'))
      p = VlnPlot(tiss, genes, group.by = group.by, do.return=TRUE, size.x.use=8, size.y.use=8, nCol=3)
      ggsave(filename, width = 8, height = 11)
      dev.off()
      }


    }
}





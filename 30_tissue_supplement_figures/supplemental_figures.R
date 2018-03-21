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
library(scales)
library(gridExtra)
library(lemon)
library(RColorBrewer)

# Commonly used variables by all Rmd files
standard.group.bys = c("cell_ontology_class", "free_annotation", 'cluster.ids')
prefix = 'allcells'

# Only use integers for xticklabels to save horizontal space
integer_breaks = function(x)
  unique(floor(pretty(seq(0, (
    max(x) + 1
  ) * 1.1))))

CHUNKSIZE = 20


# Dotplot palette
palette = brewer.pal(3, 'YlGnBu')
dotplot.cols.use = c(palette[1], palette[3])

# # Extract legend from a ggplot
# # Stolen from https://gist.github.com/crsh/be88be19233f1df4542aca900501f0fb
# # Who stole it from: http://stackoverflow.com/a/12041779/914024
# gglegend <- function(x){
#   tmp <- ggplot_gtable(ggplot_build(x))
#   leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
#   tmp$grobs[[leg]]
# }

make_filename = function(save_folder, prefix, group.by, plottype, format='pdf'){
  # Replace dots with dashes so filesystems are happy
  suffix = sub('.', '-', group.by, fixed=TRUE)
  filename = file.path(save_folder,
                       paste(prefix, suffix,
                             paste0(plottype, '.', format), sep = '_'))
  return(filename)
}

# Make all supplemental figures for a tiss object
dot_tsne_violin = function(tiss,
                           genes_to_check,
                           save_folder,
                           prefix,
                           group.bys,
                           method = 'facs') {
  if (method == 'facs') {
    expression_unit = 'log(CPM)'
  } else {
    expression_unit = 'log(CP10k)'
  }
  
  for (group.by in group.bys) {
    write(paste('group.by:', group.by, '   prefix:', prefix), stderr())
    
    # If this column is all NAs, then skip it
    if (all(is.na(tiss@meta.data[, group.by]))) {
      next
    }
    
    # Get number of colors to use
    if (group.by == 'cell_ontology_class') {
      n_annotations = dim(unique(tiss@meta.data[group.by]))[1]
      if (n_annotations > 8) {
        colors.use = c(brewer.pal(8, 'Set2'), brewer.pal(max(n_annotations - 8, 3), 'Dark2'))
      } else{
        colors.use = brewer.pal(max(n_annotations, 3), 'Set2')
      }
    } else {
      colors.use = NULL
    }
    
    # TSNE only
    filename = make_filename(save_folder, prefix, group.by, 'tsneplot')
    p = TSNEPlot(
      object = tiss,
      do.return = TRUE,
      group.by = group.by,
      no.axes = TRUE,
      pt.size = 1,
      no.legend = TRUE,
      colors.use = colors.use
    )
    ggsave(filename, width = 2, height = 2)
    # dev.off()
    
    # Legend for TSNE
    filename = make_filename(save_folder, prefix, group.by, 'tsneplot_legend')
    # Plot TSNE again just to steal the legend
    p = TSNEPlot(
      object = tiss,
      do.return = TRUE,
      group.by = group.by,
      no.axes = TRUE,
      pt.size = 1,
      no.legend = FALSE,
      label.size = 8,
      colors.use = colors.use
    ) + coord_fixed(ratio = 1) +
      xlab("tSNE 1") + ylab("tSNE 2")

    # Initialize an empty canvas!
    ggdraw()
    # Draw only gene 
    ggdraw(g_legend(p))
    ggsave(filename, width = 3, height = 3)
    dev.off()
    
    # Many of the strings are too long so we use numbers instead. The color
    # scheme is the same as the TSNE as long as we sort alphabetically.
    # Make a temporary copy of the column to store the original strings
    tiss@meta.data[, 'tmp'] = tiss@meta.data[, group.by]
    labels = sort(unique(tiss@meta.data[, group.by]))
    numbers = LETTERS[seq(1, length(labels))]
    tiss@meta.data[, group.by] = plyr::mapvalues(
      x = tiss@meta.data[, 'tmp'],
      from = labels,
      to = numbers
    )
    
    # Reverse the order of the identities to make it easier to read
    ident.use = rev(sort(unique(tiss@meta.data[, group.by])))
    
    # Break up the gene list into reasonable lists for legible pages
    chunked_genes = split(genes_to_check, ceiling(seq_along(genes_to_check) /
                                                    CHUNKSIZE))
    n_chunks = length(chunked_genes)
    for (name in names(chunked_genes)) {
      write(paste('\tchunk:', name), stderr())
      
      genes = chunked_genes[[name]]
      
      # Adjust canvas size so the plots are the same size no matter how many genes
      nCol = 4
      nRow = ceiling(length(genes) / nCol)
      
      # Set height of ridgeplots
      ridge_height = 4.4
      
      # Dotplot - enrichment of gene expression in group.by with dot size
        filename = make_filename(save_folder, prefix, group.by,
          paste0('dotplot_', name, '-of-', n_chunks)
      # Use rev(genes) to reverse the order because dotplot does Right to Left instead of left to right
      p = DotPlot(
        tiss,
        rev(genes),
        col.max = 2.5,
        plot.legend = T,
        do.return = T,
        group.by = group.by,
        cols.use = dotplot.cols.use,
      # ) #+ scale_y_reverse()
        # Reverse the yscale so the clusters appear in ascending instead of descending order
      ) + scale_y_discrete(limits = rev(levels(group.by)))
      ggsave(filename, width = 13.75, height = 10)
      # dev.off()
      
      # Ridgeplot - enrichment of gene expression in group.by with smoothed histograms
        filename = make_filename(save_folder, prefix, group.by,
          paste0('ridgeplot_', name, '-of-', n_chunks)
      plots = RidgePlot(
        tiss,
        genes,
        group.by = group.by,
        do.return = TRUE,
        nCol = 4,
        same.y.lims = TRUE,
        return.plotlist = TRUE,
        cols.use = colors.use
        # Reverse the yscale so the clusters appear in ascending instead of descending order
      )
      for (i in seq(1, length(plots))) {
        plots[[i]] = plots[[i]] + xlab(expression_unit) + scale_x_continuous(breaks = pretty_breaks(n = 4)) + scale_y_discrete(limits = rev(levels(group.by)))
        # + scale_y_reverse()
        #
      }
      plots.combined <- plot_grid(plotlist = plots, ncol = nCol)
      invisible(x = lapply(X = plots.combined, FUN = print))
      ggsave(filename,
             width = 16,
             height = nRow * ridge_height)
      # dev.off()
    }
    
    # Put the original values back
    tiss@meta.data[, group.by] = tiss@meta.data[, 'tmp']
  }
}

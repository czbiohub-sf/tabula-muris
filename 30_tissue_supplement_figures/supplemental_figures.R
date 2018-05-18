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

CHUNKSIZE = 12


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

add_black_na = function(p, group.by, colors.use, plottype = 'tsne') {
  if (group.by == 'cell_ontology_class') {
    if (plottype == 'tsne') {
      p = p + scale_colour_manual(values = colors.use, na.value = 'black')
    } else if (plottype == 'ridge') {
      p = p + scale_fill_manual(values = colors.use, na.value = 'black')
    }
  }
  return(p)
}


make_filename = function(save_folder,
                         prefix,
                         group.by,
                         plottype,
                         format = 'pdf') {
  # Replace dots with dashes so filesystems are happy
  suffix = sub('.', '-', group.by, fixed = TRUE)
  filename = file.path(save_folder,
                       paste(prefix, suffix,
                             paste0(plottype, '.', format), sep = '_'))
  return(filename)
}

write_caption = function(caption, pdf_name) {
  write(caption, sub('.pdf', '.txt', pdf_name))
}

get_valid_genes = function(tiss, genes_to_check) {
  genes_not_found = setdiff(genes_to_check, rownames(tiss@scale.data))
  
  if (length(genes_not_found) > 0) {
    for (gene in genes_not_found) {
      write(paste("This gene was not found in the data:", gene),
            stderr())
      similar_genes = sort(rownames(tiss@raw.data)[agrep(gene, rownames(tiss@raw.data))])
      write(paste(
        "\t",
        gene,
        "was not found, but here is a simlar name:",
        paste(similar_genes, sep = ', ')
      ),
      stderr())
      
    }
    
  }
  genes_to_check = intersect(genes_to_check, rownames(tiss@scale.data))
  return(genes_to_check)
}

violinplot_and_save = function(tiss,
                              save_folder,
                              prefix,
                              group.by,
                              genes,
                              colors.use = NULL,
                              i = 0,
                              n = 0,
                              suffix = NULL,
                              method = 'facs') {
  if (length(colors.use) == 0) {
    colors.use = get_colors(tiss, group.by)
  }
  if (method == 'facs') {
    expression_unit = 'log(CPM)'
  } else {
    expression_unit = 'log(CP10k)'
  }
  
  
  # Add letter to the string labels groupby for legends
  if (group.by != 'cluster.ids') {
    labels = sort(unique(tiss@meta.data[, group.by]))
    letters = LETTERS[seq(1, length(labels))]
    # Many of the strings are too long so we use numbers instead. The color
    # scheme is the same as the TSNE as long as we sort alphabetically.
    # Make a temporary copy of the column to store the original strings
    tiss@meta.data[, 'tmp'] = tiss@meta.data[, group.by]
    tiss@meta.data[, group.by] = plyr::mapvalues(x = tiss@meta.data[, 'tmp'],
                                                 from = labels,
                                                 to = letters)
  } else {
    # Make sure cluster ids are numeric so they're sorted properly
    tiss@meta.data[, group.by] = as.numeric(tiss@meta.data[, group.by])
  }
  
  
  if (length(suffix) == 0) {
    suffix = paste0(i, '-of-', n)
  }
  filename = make_filename(save_folder,
                           prefix,
                           group.by,
                           paste('violinplot', suffix, sep = '_'))
  # Adjust canvas size so the plots are the same size no matter how many genes
  nCol = 3
  nRow = ceiling(length(genes) / nCol)
  
  # Set height of violinplots
  ridge_height = 5.5
  
  plots = VlnPlot(
    tiss,
    genes,
    group.by = group.by,
    do.return = TRUE,
    nCol = nCol,
    same.y.lims = TRUE,
    return.plotlist = TRUE,
    cols.use = colors.use
  )
  for (i in seq(1, length(plots))) {
    p = plots[[i]]
    p = p + coord_flip()
      p = p + xlab(expression_unit) + scale_x_continuous(breaks = pretty_breaks(n = 4)) + scale_y_discrete(limits = rev(levels(group.by)))
    p = add_black_na(p, group.by, colors.use, plottype = 'ridge')
    plots[[i]] = p
  }
  plots.combined <- plot_grid(plotlist = plots, ncol = nCol)
  invisible(x = lapply(X = plots.combined, FUN = print))
  ggsave(filename,
         width = 16,
         height = nRow * ridge_height)
  
  genes_filename = sub('.pdf', '_genes.txt', filename)
  write(c(genes[1], genes[length(genes)]), genes_filename)
  
  
  return(filename)
}


tsne_and_legend = function(tiss,
                           save_folder,
                           prefix,
                           group.by,
                           colors.use = NULL,
                           do.label = FALSE) {
  if (length(colors.use) == 0) {
    colors.use = get_colors(tiss, group.by)
  }
  # TSNE only
  tsne_filename = make_filename(save_folder, prefix, group.by, 'tsneplot')
  p = TSNEPlot(
    object = tiss,
    do.return = TRUE,
    group.by = group.by,
    no.axes = TRUE,
    pt.size = 1,
    no.legend = TRUE,
    colors.use = colors.use,
    do.label = do.label
  )
  p = add_black_na(p, group.by, colors.use, plottype = 'tsne')
  ggsave(tsne_filename, width = 4, height = 4)
  # dev.off()
  
  # Add letter to the string labels groupby for legends
  if (group.by != 'cluster.ids') {
    labels = sort(unique(tiss@meta.data[, group.by]))
    letters = LETTERS[seq(1, length(labels))]
    new_labels = paste0('(', letters, ') ', labels)
    group.by.plus = paste0(group.by, '_letter')
    tiss@meta.data[, group.by.plus] = plyr::mapvalues(x = tiss@meta.data[, group.by],
                                                      from = labels,
                                                      to = new_labels)
  } else {
    # Make sure cluster ids are numeric so they're sorted properly
    tiss@meta.data[, group.by] = as.numeric(tiss@meta.data[, group.by])
    group.by.plus = group.by
  }
  
  
  legend_filename = make_filename(save_folder, prefix, group.by, 'tsneplot_legend')
  # Plot TSNE again just to steal the legend
  p = TSNEPlot(
    object = tiss,
    do.return = TRUE,
    group.by = group.by.plus,
    no.axes = TRUE,
    pt.size = 1,
    no.legend = FALSE,
    label.size = 8,
    colors.use = colors.use
  ) + coord_fixed(ratio = 1) +
    xlab("tSNE 1") + ylab("tSNE 2")
  p = add_black_na(p, group.by, colors.use, plottype = 'tsne')
  
  # Initialize an empty canvas!
  ggdraw()
  # Draw only the legend
  ggdraw(g_legend(p))
  ggsave(legend_filename, width = 8, height = 4)
  dev.off()
  
  return(tsne_filename)
}

get_colors = function(tiss, group.by) {
  annotations = sort(unique(tiss@meta.data[, group.by]))
  n_annotations = length(annotations)
  
  # Get number of colors to use
  if (group.by == 'cell_ontology_class') {
    if (n_annotations > 8) {
      if (n_annotations > 16) {
        colors.use = c(brewer.pal(8, 'Set2'),
                       brewer.pal(8, 'Dark2'),
                       brewer.pal(max(n_annotations - 16, 3), 'Pastel2'))
        
      } else {
        colors.use = c(brewer.pal(8, 'Set2'), brewer.pal(max(n_annotations - 8, 3), 'Dark2'))
      }
    } else{
      colors.use = brewer.pal(max(n_annotations, 3), 'Set2')
    }
    # Add black for NA
    # colors.use = c(colors.use, 'black')
  } else {
    colors.use = NULL
  }
  return (colors.use)
}

dotplot_and_save = function(tiss,
                            save_folder,
                            prefix,
                            group.by,
                            genes,
                            i = 0,
                            n = 0,
                            suffix = NULL,
                            method = 'facs') {
  # Add letter to the string labels groupby for legends
  if (group.by != 'cluster.ids') {
    labels = sort(unique(tiss@meta.data[, group.by]))
    letters = LETTERS[seq(1, length(labels))]
    # Many of the strings are too long so we use numbers instead. The color
    # scheme is the same as the TSNE as long as we sort alphabetically.
    # Make a temporary copy of the column to store the original strings
    tiss@meta.data[, 'tmp'] = tiss@meta.data[, group.by]
    tiss@meta.data[, group.by] = plyr::mapvalues(x = tiss@meta.data[, 'tmp'],
                                                 from = labels,
                                                 to = letters)
  } else {
    # Make sure cluster ids are numeric so they're sorted properly
    tiss@meta.data[, group.by] = as.numeric(tiss@meta.data[, group.by])
  }
  
  
  if (length(suffix) == 0) {
    suffix = paste0(i, '-of-', n)
  }
  
  # Dotplot - enrichment of gene expression in group.by with dot size
  filename = make_filename(save_folder,
                           prefix,
                           group.by,
                           paste('dotplot', suffix, sep = '_'))
  # Use rev(genes) to reverse the order because dotplot does Right to Left instead of left to right
  p = DotPlot(
    tiss,
    rev(genes),
    col.max = 2.5,
    plot.legend = T,
    do.return = T,
    group.by = group.by,
    cols.use = dotplot.cols.use,
  )
  # ) + scale_y_discrete(trans='reverse')
  # ) #+ scale_y_reverse()
  #   # Reverse the yscale so the clusters appear in ascending instead of descending order
  # ) + scale_y_discrete(limits = rev(levels(ident.use)))
  ggsave(filename, width = 13.75, height = 10)
  genes_filename = sub('.pdf', '_genes.txt', filename)
  write(c(genes[1], genes[length(genes)]), genes_filename)
  return(filename)
}

# Make all supplemental figures for a tiss object
dot_tsne_ridge = function(tiss,
                          genes_to_check,
                          save_folder,
                          prefix,
                          group.bys,
                          method = 'facs') {
  genes_to_check = get_valid_genes(tiss, genes_to_check)
  
  for (group.by in group.bys) {
    write(paste('group.by:', group.by, '   prefix:', prefix), stderr())
    
    # If this column is all NAs, then skip it
    if (all(is.na(tiss@meta.data[, group.by])) ||
        all(tiss@meta.data[, group.by] == "NA")) {
      next
    }
    
    annotations = sort(unique(tiss@meta.data[, group.by]))
    n_annotations = length(annotations)
    
    write(paste('Annotations in', tiss@project.name, method, group.by, ':'),
          stderr())
    write(paste0('\t', annotations), stderr())
    
    # Get number of colors to use
    if (group.by == 'cell_ontology_class') {
      if (n_annotations > 8) {
        if (n_annotations > 16) {
          colors.use = c(
            brewer.pal(8, 'Set2'),
            brewer.pal(8, 'Dark2'),
            brewer.pal(max(n_annotations - 16, 3), 'Pastel2')
          )
          
        } else {
          colors.use = c(brewer.pal(8, 'Set2'), brewer.pal(max(n_annotations - 8, 3), 'Dark2'))
        }
      } else{
        colors.use = brewer.pal(max(n_annotations, 3), 'Set2')
      }
      # Add black for NA
      # colors.use = c(colors.use, 'black')
    } else {
      colors.use = NULL
    }
    
    tsne_and_legend(tiss, save_folder, prefix, group.by, colors.use)
    
    
    
    
    # Break up the gene list into reasonable lists for legible pages
    chunked_genes = split(genes_to_check, ceiling(seq_along(genes_to_check) /
                                                    CHUNKSIZE))
    n_chunks = length(chunked_genes)
    for (i in names(chunked_genes)) {
      write(paste('\tchunk:', i), stderr())
      
      genes = chunked_genes[[i]]
      
      dotplot_and_save(
        tiss,
        save_folder,
        prefix,
        group.by,
        genes,
        i = i,
        n = n_chunks,
        method = method
      )
      
      violinplot_and_save(
        tiss,
        save_folder,
        prefix,
        group.by,
        genes,
        colors.use,
        i = i,
        n = n_chunks,
        method = method
      )
      # dev.off()
    }
    
    
  }
}

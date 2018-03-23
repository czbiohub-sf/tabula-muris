#!/usr/bin/env python3.6
# coding: utf-8

import glob
import os
import string

import click
import pandas as pd

FIGURE_ORDER = 'tsneplot', 'ridgeplot', 'dotplot'
GROUPBY_ORDER = 'cell_ontology_class', 'cluster-ids', 'free_annotation'
FIGURE_FOLDER = '30_tissue_supplement_figures'
PATTERN = '^(?P<subset>[a-zA-Z\d]+)_(?P<groupby>[\w\-]+)_(?P<plottype>[a-z]+plot)(_?:(?P<i>\d+)\-of\-(?P<n>\d+))?_?(?P<is_legend>legend)?.pdf$'


class FigureTeX:

    def __init__(self, plottype, tissue, method, subset, groupby, i=None,
                 n=None, labels=None):
        self.plottype = plottype
        self.tissue = tissue
        self.method = method
        self.subset = subset
        self.groupby = groupby
        self.i = i
        self.n = n
        self.labels = labels

    @property
    def is_iterative(self):
        return self.i is not None and self.i > 1

    @property
    def plottype_tex(self):
        if self.plottype == 'tsneplot':
            return 't-Distributed stochastic neighbor embedding (tSNE) plot'
        else:
            # ridgeplot and dotplot
            return self.plottype.capitalize()

    @property
    def plottype_title(self):
        if self.plottype == 'tsneplot':
            return 't-SNE plot'
        else:
            # ridgeplot and dotplot
            return self.plottype.capitalize()

    @property
    def subset_tex(self):
        if self.subset == 'allcells':
            return "all cells"
        else:
            return self.subset

    @property
    def groupby_tex(self):
        """TeX-formatted groupby"""
        return self.groupby.replace('_', ' ').replace('.', ' ')

    @property
    def plot_shows(self):
        if self.plottype == 'tsneplot':
            return ''
        else:
            # Dotplot and ridgeplot
            return ' showing gene expression enrichment'

    @property
    def caption_start(self):
        if self.plottype == 'tsneplot':
            return "Top,"
        else:
            # ridgeplot and dotplot
            return ''

    @property
    def caption_end(self):
        if self.plottype == 'tsneplot':
            groupby_tex = self.groupby_tex
            if self.groupby != 'cluster-ids':
                groupby_tex += ' (and letter abbreviation)'

            return f"Bottom, legend mapping {groupby_tex} to colors"
        else:
            # ridgeplot and dotplot
            letter_labels = ', '.join([f'{letter}: {label}' for letter, label
                                       in zip(string.ascii_uppercase,
                                              self.labels)])
            letter_labels += '.'
            return letter_labels

    @property
    def caption(self):
        words = [self.caption_start, self.plottype_tex]
        if self.is_iterative:
            words.append(f'({self.i} of {self.n})')
        words.extend([self.plot_shows, 'in', self.groupby_tex, 'of',
                      self.tissue, self.method + '.', self.caption_end])
        sentence = ' '.join(words)
        return f'\caption{{{sentence}}}'

    @property
    def graphics_options(self):
        if self.plottype == 'tsneplot':
            return 'height=.35\\textheight'
        if self.plottype == 'ridgeplot':
            return 'width=.75\\textwidth'
        if self.plottype == 'dotplot':
            return 'angle=90, height=.7\\textheight'

    @property
    def pdf(self):
        prefix = f'{self.subset}_{self.groupby}_{self.plottype}'
        if self.is_iterative:
            basename = f'{prefix}_{self.i}-of-{self.n}.pdf'
        else:
            basename = 'f{prefix}.pdf'
        return os.path.join('..', '30_tissue_supplement_figures', self.tissue,
                            self.method, basename)

    @property
    def legend(self):
        if self.plottype == 'tsneplot':
            pdf = self.pdf.replace('.pdf', '_legend.pdf')
            return f'\includegraphics[{self.graphics_options}]{{{pdf}}}'
        else:
            return ''

    @property
    def section_title(self):
        title = self.plottype_title
        if self.is_iterative:
            title += f' {self.i} of {self.n}}}'
        return title

    @property
    def figure_code(self):
        code = f"""\\newpage
\subsection{{{self.section_title}}}
\\begin{{figure}}[h]
\centering
\includegraphics[{self.graphics_options}]{{{self.pdf}}}{self.legend}
{self.caption}
\end{{figure}}
"""
        return code


@click.command()
@click.option('--tissue', default='all')
@click.option('--method', default='all')
def cli(tissue, method):
    tissue = '*' if tissue == 'all' else tissue.capitalize()
    method = '*' if method == 'all' else method.lower()
    globber = os.path.join('..', FIGURE_FOLDER, tissue, method)
    tissue_method_paths = glob.glob(globber)
    for tissue_method_path in tissue_method_paths:
        figures = glob.glob(os.path.join(tissue_method_path, '*.pdf'))
        if len(figures) == 0:
            continue
        tissue_path, method = os.path.split(tissue_method_path)
        figure_path, tissue = os.path.split(tissue_path)
        print(f'\n--- tissue: "{tissue}", method: "{method}" ---')
        basename = '_'.join([tissue, method, 'annotation.csv'])
        annotation = pd.read_csv(os.path.join('..', '00_data_ingest',
                                              '03_tissue_annotation_csv',
                                              basename))

        basenames = [os.path.basename(f) for f in figures]
        basenames = pd.Series(basenames, index=basenames, name='basename')
        parameters = basenames.str.extractall(PATTERN)
        print(parameters)



if __name__ == "__main__":
    cli()

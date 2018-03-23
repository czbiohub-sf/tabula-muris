#!/usr/bin/env python3.6
# coding: utf-8

import glob
import os
import string

import click
import pandas as pd

SUBSET_ORDER = 'allcells',
PLOT_ORDER = 'tsneplot', 'ridgeplot', 'dotplot'
GROUPBY_ORDER = 'cell_ontology_class', 'cluster-ids', 'free_annotation'
ORDER_DEFAULTS = {'subset': SUBSET_ORDER, 'plottype': PLOT_ORDER,
                  'groupby': GROUPBY_ORDER}

FIGURE_FOLDER = '30_tissue_supplement_figures'
PATTERN = '^(?P<subset>[a-zA-Z\d]+)_(?P<groupby>[\w\-]+)_(?P<plottype>[a-z]+plot)(_?:(?P<i>\d+)\-of\-(?P<n>\d+))?_?(?P<extra>[a-z\-A-Z0-9]+)?.pdf$'
FRONTMATTER = r"""\documentclass[11pt,]{article}   	% use "amsart" instead of "article" for AMSLaTeX format
\usepackage[a4paper, margin=1.75cm]{geometry}
\usepackage[page]{totalcount}
%\geometry{landscape}                		% Activate for rotated page geometry
%\usepackage[parfill]{parskip}    		% Activate to begin paragraphs with an empty line rather than an indent
\usepackage{graphicx}				% Use pdf, png, jpg, or epsÂ§ with pdflatex; use eps in DVI mode
								% TeX will automatically convert eps --> pdf in pdflatex		
\usepackage{amssymb}

% Reference sections by name
\usepackage{nameref}

% Make command to reference section name
\makeatletter
\newcommand*{\currentsection}{\@currentlabelname}
\makeatother

\renewcommand{\sectionmark}[1]{\markboth{#1}{}} % set the \leftmark


% Rotated figures
\usepackage{rotating}

% Add header
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\rhead{TISSUE METHOD Figure Packet}
\lhead{\itshape\nouppercase{\leftmark}}
\rfoot{Page \thepage~of~\totalpages}

%SetFonts

%SetFonts


\title{TISSUE METHOD Figure Packet}
%\author{The Author}
\date{}							% Activate to display a given date or no date

\begin{document}
\maketitle
\tableofcontents
%\listoffigures
"""

SECTION = r"""\newpage
\section{SUBSET, labeled by GROUPBY}
"""

ENDMATTER = r"\end{document}"


class FigureTeXGenerator:
    def __init__(self, pdf, plottype, tissue, method, subset, groupby, i=None,
                 n=None, labels=None):
        self.pdf = pdf
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

    def generate_code(self):
        code = f"""
\\newpage
\subsection{{{self.section_title}}}
\\begin{{figure}}[h]
\centering
\includegraphics[{self.graphics_options}]{{{self.pdf}}}{self.legend}
{self.caption}
\end{{figure}}

"""
        return code


def get_category_order(column, defaults):
    column_unique = set(column.astype(str).unique())
    print(column_unique)
    remaining = column_unique.difference(defaults)
    categories = list(defaults) + list(remaining)
    return categories


def add_categorical_order(parameters, cols=('subset', 'groupby', 'plottype')):
    """Ensures first group is always TSNE of allcells + cell_ontology_class"""
    for col in cols:
        defaults = ORDER_DEFAULTS[col]
        categories = get_category_order(parameters[col], defaults)
        parameters[col] = pd.Categorical(parameters[col],
                                         categories=categories)
    return parameters


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
        tex = FRONTMATTER.replace('TISSUE', tissue).replace("METHOD", method)
        print(f'\n--- tissue: "{tissue}", method: "{method}" ---')
        basename = '_'.join([tissue, method, 'annotation.csv'])
        annotation = pd.read_csv(os.path.join('..', '00_data_ingest',
                                              '03_tissue_annotation_csv',
                                              basename))

        basenames = [os.path.basename(f) for f in figures]
        basenames = pd.Series(basenames, index=basenames, name='basename')
        parameters = basenames.str.extractall(PATTERN)
        parameters.index = parameters.index.droplevel(-1)
        parameters = add_categorical_order(parameters)
        # print(parameters)

        # Remove legend figures because they're auto-added
        ind = parameters['extra'] != 'legend'
        parameters = parameters.query('extra != "legend"')
        grouped = parameters.groupby(['subset', 'groupby', 'plottype'])
        for (subset, groupby, plottype), row in grouped:

            # Replaced dots with dashes for filenames, need to change back
            # for column referencing
            groupby_col = groupby.replace('-', '.')
            print(row)
            groupby_unique = annotation[groupby_col].astype(str).unique()
            print(groupby_unique)
            labels = sorted(groupby_unique)
            tex += SECTION.replace('SUBSET', subset).replace("GROUPBY",
                                                             groupby)
            # This groupby iterates by row but we still need to grab the first
            # item because pandas doesn't cast the row to a vector
            pdf = row.index[0]
            i = row['i'].iloc[0] if pd.notnull(row['i']).all() else None
            n = row['n'].iloc[0] if pd.notnull(row['n']).all() else None
            figuretexgen = FigureTeXGenerator(pdf, plottype, tissue,
                                              method, subset, groupby,
                                              i=i, n=n,
                                              labels=labels)
            figuretex = figuretexgen.generate_code()
            tex += figuretex

        tex += ENDMATTER
        filename = f'{tissue}_{method}_auto_generated.tex'
        with open(filename, 'w') as f:
            f.write(tex)


if __name__ == "__main__":
    cli()

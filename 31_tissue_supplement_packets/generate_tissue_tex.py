#!/usr/bin/env python3.6
# coding: utf-8

import glob
import os
import string

import click
import pandas as pd
import yaml


SUBSET_ORDER = 'allcells',
PLOT_ORDER = 'tsneplot', 'ridgeplot', 'dotplot'
GROUPBY_ORDER = 'cell_ontology_class', 'cluster-ids', 'free_annotation'
ORDER_DEFAULTS = {'subset': SUBSET_ORDER, 'plottype': PLOT_ORDER,
                  'groupby': GROUPBY_ORDER}

FIGURE_FOLDER = '30_tissue_supplement_figures'
PATTERN = '^(?P<subset>[a-zA-Z\d]+)_(?P<groupby>[\w\-]+)_(?P<plottype>[a-z]+plot)(_?:(?P<i>\d+)\-of\-(?P<n>\d+))?_?(?P<extra>[a-z\-A-Z0-9]+)?.pdf$'

SUBSECTION = r"""
\newpage
\subsection{SUBSET, labeled by GROUPBY}
"""

ENDMATTER = r"\end{document}"



def method_tex(method):
    if method == 'droplet':
        return method.capitalize()
    else:
        return method.replace('_', ' ').upper()


class TeXGenerator:
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
    def tissue_tex(self):
        return self.tissue.replace('_', ' ')

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
    def method_tex(self):
        return method_tex(self.method)

    @property
    def subset_tex(self):
        if self.subset.lower().endswith('cells'):
            subset = self.subset.lower().split('cells')[0]
            return subset + ' cells'
        else:
            return self.subset

    @property
    def groupby_tex(self):
        """TeX-formatted groupby"""
        return self.groupby.replace('_', ' ').replace('.', ' ').replace('-', ' ')

    @property
    def subsection_tex(self):
        tex = SUBSECTION.replace('GROUPBY', self.groupby_tex.title()).replace(
            "SUBSET", self.subset_tex.title())
        return tex

    @property
    def labels_tex(self):
        tex = [str(x).replace('_', ' ') for x in self.labels]
        return tex


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
                                              self.labels_tex)])
            letter_labels += '.'
            return letter_labels

    @property
    def caption(self):
        caption_file = self.pdf.replace('.pdf', '.txt')
        if os.path.exists(caption_file):
            with open(caption_file) as f:
                sentence = f.read()
        else:
            words = [self.caption_start, self.plottype_tex]
            if self.is_iterative:
                words.append(f'({self.i} of {self.n})')
            words.extend([self.plot_shows, 'in', self.groupby_tex, 'of',
                          self.tissue_tex, self.method_tex + '.', self.caption_end])
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
        if self.plottype == 'tsneplot' and not('expression' in self.groupby):
            pdf = self.pdf.replace('.pdf', '_legend.pdf')
            return f'\includegraphics[{self.graphics_options}]{{{pdf}}}'
        else:
            return ''

    @property
    def subsubsection_title(self):
        title = self.plottype_title
        if self.is_iterative:
            title += f' {self.i} of {self.n}}}'
        return title

    @property
    def figure_tex(self):
        code = f"""
\subsubsection{{{self.subsubsection_title}}}
\\begin{{figure}}[h]
\centering
\includegraphics[{self.graphics_options}]{{{self.pdf}}}
{self.legend}
{self.caption}
\end{{figure}}

"""
        return code

    def make_table_tex(self, counts):
        tex = f'\subsubsection{{Table of cell counts in {self.subset_tex}, per {self.groupby_tex}}}'
        tex += r"""\begin{table}[h]
\centering
\label{my-label}
\begin{tabular}{@{}ll@{}}
\toprule
"""
        tex += f'\n{self.groupby_tex.capitalize()}& Number of cells \\\\ \midrule'
        for label, n in zip(counts.index, counts.values):
            try:
                label = label.replace("_", " ")
            except AttributeError:
                # Label is an int
                pass
            tex += f'\n{label} & {n} \\\\\n'
        tex += r'\bottomrule'
        tex += f'''
\end{{tabular}}
\caption{{Cell counts for {self.subset_tex}, per {self.groupby_tex}.}}
\end{{table}}
'''
        return tex


def get_category_order(column, defaults):
    column_unique = set(column.astype(str).unique())
    # print(column_unique)
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
    tissue_method_paths = sorted(glob.glob(globber))

    tex_filenames = []
    for tissue_method_path in tissue_method_paths:
        figures = glob.glob(os.path.join(tissue_method_path, '*.pdf'))
        if len(figures) == 0:
            continue
        tissue_path, method = os.path.split(tissue_method_path)
        figure_path, tissue = os.path.split(tissue_path)
        tissue_tex = tissue.replace('_', ' ')

        basename_yaml = f'{tissue.lower()}_{method}.yaml'
        filename_yaml = os.path.join('..', '28_tissue_yamls_for_supplement', basename_yaml)
        with open(filename_yaml) as f:
            yaml_data = yaml.load(f)

        tex = f'''\\newpage
\section{{{tissue_tex} {method_tex(method)}}}
'''

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
        parameters = parameters.query('extra != "legend"')
        grouped = parameters.groupby(['subset', 'groupby', 'plottype'])

        # Section counter
        j = 0
        prev_subset = ''
        prev_groupby = ''
        for (subset, groupby, plottype), row in grouped:
            # Replaced dots with dashes for filenames, need to change back
            # for column referencing
            groupby_col = groupby.replace('-', '.')
            try:
                try:
                    groupby_unique = annotation[groupby_col].unique()
                    labels = sorted(groupby_unique)
                except TypeError:
                    groupby_unique = annotation[groupby_col].astype(str).unique()
                    labels = sorted(groupby_unique)
            except KeyError:
                labels = None
            # This groupby iterates by row but we still need to grab the first
            # item because pandas doesn't cast the row to a vector
            pdf = os.path.join(tissue_method_path, row.index[0])
            i = row['i'].iloc[0] if pd.notnull(row['i']).all() else None
            n = row['n'].iloc[0] if pd.notnull(row['n']).all() else None
            tex_generator = TeXGenerator(pdf, plottype, tissue,
                                        method, subset, groupby,
                                        i=i, n=n,
                                        labels=labels)

            # Weird hack for adding sections
            if prev_groupby != groupby or prev_subset != subset:
                j = 0
            else:
                j += 1

            print(f'\tsubset: {subset}, groupby: {groupby}, plottype: {plottype}, j: {j}, i: {i}, n: {n}')

            if j == 0:
                tex += tex_generator.subsection_tex

                if subset != 'allcells':
                    try:
                        subset_yaml = yaml_data['SUBSET'][groupby.upper()]
                        print(subset_yaml)
                        # Add table of counts
                        column = subset_yaml['FILTER_COLUMN']
                        value = subset_yaml['FILTER_VALUE']
                    except KeyError:
                        column = 'cell_ontology_class'
                        value = subset.lower().replace('_', ' ').rstrip('s')
                    subset_annotation = annotation.query(f'{column} == "{value}"')
                else:
                    subset_annotation = annotation

                if 'expression' not in groupby_col:
                    counts = subset_annotation.groupby(groupby_col).size()
                    tex += tex_generator.make_table_tex(counts)
                    tex += r'''
\newpage'''

            else:
                tex += r'''
\newpage'''

            # Add section title and figure tex
            tex += tex_generator.figure_tex

            prev_subset = subset
            prev_groupby = groupby

        filename = f'{tissue}_{method}_auto_generated.tex'
        with open(filename, 'w') as f:
            f.write(tex)
        tex_filenames.append(filename)

    print(tex_filenames)
    with open('tissue_supplement_template.tex') as f:
        content = f.read()

    include_commands = '\n'.join(['\include{' + x.split('.tex')[0] + '}'
                                  for x in tex_filenames])
    content = content.replace('%%% --- INCLUDE TISSUES HERE --- %',
                              include_commands)

    with open('tissue_supplement.tex', 'w') as g:
        g.write(content)



if __name__ == "__main__":
    cli()

#!/usr/bin/env python3.6
# coding: utf-8

import glob
from functools import cmp_to_key
import locale
import os
import re
import string
import sys

import click
import pandas as pd
import yaml



INT_RE = re.compile(r"^[-]?\d+$")

def is_int(s):
    """String contains integer"""
    return INT_RE.match(str(s)) is not None

SUBSET_ORDER = 'allcells',
PLOT_ORDER = 'tsneplot', 'ridgeplot', 'dotplot'
GROUPBY_ORDER = ('highlighted', 'cell_ontology_class', 'cluster-ids',
                 'free_annotation')
ORDER_DEFAULTS = {'subset': SUBSET_ORDER, 'plottype': PLOT_ORDER,
                  'groupby': GROUPBY_ORDER}

FIGURE_FOLDER = '30_tissue_supplement_figures'
PATTERN = '^(?P<subset>[a-zA-Z\d]+)_(?P<groupby>[\w\-]+)_(?P<plottype>[a-z]+plot)(_?:(?P<i>\d+)\-of\-(?P<n>\d+))?_?(?P<extra>[a-z\-A-Z0-9]+)?.pdf$'

SUBSECTION = r"""
\subsection{SUBSET, labeled by GROUPBY}
"""

ENDMATTER = r"\end{document}"


def method_tex(method):
    if method == 'droplet':
        return method.capitalize()
    else:
        return method.replace('_', ' ').upper()


def alphabetical_sort(values):
    if all(isinstance(x, str) for x in values):
        return sorted(values, key=lambda x: x.lower())
    else:
        return sorted(values)


def unique_sorted(series):
    if series.notnull().all():
        groupby_unique = series.unique()
        labels = alphabetical_sort(groupby_unique)
    else:
        # There's NA and we need to add it separately
        groupby_unique = series.dropna().unique()
        labels = alphabetical_sort(groupby_unique) + ['NA']
    return labels


class TeXGenerator:
    def __init__(self, pdf, plottype, tissue, method, subset, groupby, i=None,
                 n=None, labels=None, extra=None):
        self.pdf = pdf
        self.plottype = plottype
        self.tissue = tissue
        self.method = method
        self.subset = subset
        self.groupby = groupby
        self.i = i
        self.n = n
        self.labels = labels
        self.extra = extra

    @property
    def is_iterative(self):
        return self.n is not None and self.n > 1

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
        elif 'subset' in self.subset.lower():
            split = self.subset.lower().split('subset')
            subset = 'Subset ' + split[-1].upper()
            return subset
        else:
            return self.subset

    @property
    def groupby_tex(self):
        """TeX-formatted groupby"""
        groupby = self.groupby.replace('_', ' ').replace('.', ' ').\
            replace('-', ' ')
        if 'Subset' in groupby:
            groupby = 'Subset' + groupby.split('Subset')[-1].upper()
        if groupby.lower() != 'cluster ids':
            groupby = groupby.title()
        else:
            groupby = 'Cluster IDs'
        return f'\emph{{{groupby}}}'

    @property
    def subsection_tex(self):
        if self.groupby != 'highlighted':
            tex = SUBSECTION.replace('GROUPBY', self.groupby_tex)
            tex = tex.replace("SUBSET", self.subset_tex.title())
        else:
            tex = f'\subsection{{{self.subset_tex.title()}, ' \
                  f'highlighted from All Cells tSNE}}'
        return tex

    @property
    def labels_tex(self):
        tex = [str(x).replace('_', ' ') for x in self.labels]
        return tex

    @property
    def genes_tex(self):
        """First and last gene in plot, if applicable"""
        filename = self.pdf.replace('.pdf', '_genes.txt')
        if os.path.exists(filename):
            with open(filename) as f:
                lines = [x.strip() for x in f.readlines()]
            first_gene = lines[0]
            last_gene = lines[-1]
            tex = f', {first_gene}-{last_gene}'
            return tex
        else:
            return ''

    @property
    def labels_are_ints(self):
        """Boolean if labels can be cast as int"""
        return all(map(is_int, self.labels))


    @property
    def letter_to_label(self):
        """For ridgeplot and dotplot when groupby is not cluster ids (ints)"""
        if not self.labels_are_ints:
            # ridgeplot and dotplot
            letter_labels = ', '.join([f'{letter}: {label}' if
                                       label != 'NA' else f'{label}: Unknown'
                                       for letter, label
                                       in zip(string.ascii_uppercase,
                                              self.labels_tex)])
            letter_labels += '.'
            return letter_labels
        else:
            return ''

    @property
    def extra_tex(self):
        tex = self.extra.replace('-', ' ').replace('_', ' ').title()
        return tex


    @property
    def plot_shows(self):
        if self.plottype == 'tsneplot':
            return ''
        else:
            # Dotplot and ridgeplot
            return ' showing gene expression enrichment in'

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
            return self.letter_to_label

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
            words.extend([self.plot_shows,
                          self.groupby_tex, 'labels in', self.subset_tex,
                          'of',
                          self.tissue_tex, self.method_tex + '.',
                          self.caption_end])
            sentence = ' '.join(words)
        return f'\caption{{{sentence}}}'

    @property
    def graphics_options(self):
        if self.plottype == 'tsneplot':
            return 'height=.35\\textheight'
        if self.plottype == 'ridgeplot':
            return 'width=.7\\textwidth'
        if self.plottype == 'dotplot':
            return 'angle=90, height=.6\\textheight'
        else:
            return 'width=.6\\textwidth'

    @property
    def legend(self):
        if self.plottype == 'tsneplot' and not('expression' in self.groupby):
            pdf = self.pdf.replace('.pdf', '_legend.pdf')
            if os.path.exists(pdf):
                return f'\includegraphics[{self.graphics_options}]{{{pdf}}}'
            else:
                return ''
        else:
            return ''

    @property
    def subsubsection_title(self):
        title = self.plottype_title
        if self.is_iterative:
            title += f' ({self.i} of {self.n}{self.genes_tex})'
        elif self.extra is not None and self.extra:
            title += f' ({self.extra_tex})'
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
        tex += f'\n{self.groupby_tex}& Number of cells \\\\ \midrule'
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
@click.option('--figure-folder', default=FIGURE_FOLDER)
@click.option('--tissue', default='all')
@click.option('--method', default='all')
def cli(figure_folder, tissue, method):
    """Create per-tissue tex files
    
    Total number of tissue figures: 770
    """
    tissue = '*' if tissue == 'all' else tissue
    method = '*' if method == 'all' else method.lower()
    globber = os.path.join('..', figure_folder, tissue, method)
    print(f'globber: {globber}')
    tissue_method_paths = sorted(glob.glob(globber))

    tex_filenames = []
    for tissue_method_path in tissue_method_paths:
        print(f'\n\ntissue_method_path: {tissue_method_path}')
        figures = glob.glob(os.path.join(tissue_method_path, '*.pdf'))
        figures_str = "\n\t\t".join(figures)
        print(f'\tFigures: {figures_str}')
        if len(figures) == 0:
            continue
        tissue_path, method = os.path.split(tissue_method_path)
        figure_path, tissue = os.path.split(tissue_path)
        tissue_tex = tissue.replace('_', ' ')

        basename_yaml = f'{tissue.lower()}_{method}.yaml'
        filename_yaml = os.path.join('..', '28_tissue_yamls_for_supplement',
                                     basename_yaml)
        try:
            with open(filename_yaml) as f:
                yaml_data = yaml.load(f)
        except FileNotFoundError:
            # Microbiome doesn't have a yaml
            yaml_data = None

        tex = f'''\\clearpage
\section{{{tissue_tex} {method_tex(method)}}}
'''

        print(f'\n--- tissue: "{tissue}", method: "{method}" ---')
        if tissue != 'Microbiome':
            basename = '_'.join([tissue, method, 'annotation.csv'])
            annotation = pd.read_csv(os.path.join('..', '00_data_ingest',
                                                  '03_tissue_annotation_csv',
                                                  basename))
        else:
            # Microbiome doesn't have a per-cell annotation
            annotation = None

        basenames = [os.path.basename(f) for f in figures]
        basenames = pd.Series(basenames, index=basenames, name='basename')
        parameters = basenames.str.extractall(PATTERN)
        parameters.index = parameters.index.droplevel(-1)
        parameters = add_categorical_order(parameters)
        # print(parameters)

        # Remove legend figures because they're auto-added
        parameters = parameters.query('extra != "legend"')
        parameters.loc[:, 'extra'] = parameters['extra'].fillna('')
        grouped = parameters.groupby(['subset', 'groupby', 'plottype', 'extra'])

        # Section counterr
        j = 0
        prev_subset = ''
        prev_groupby = ''
        figures_added = []

        for k, ((subset, groupby, plottype, extra), row) in enumerate(grouped):
            # if extra == 'allcells':
            #     continue

            # Replaced dots with dashes for filenames, need to change back
            # for column referencing
            if k > 0 :
                tex += r'''
\clearpage
'''

            groupby_col = groupby.replace('-', '.')
            try:
                labels = unique_sorted(annotation[groupby_col])
            except (KeyError, TypeError):
                # KeyError: groupby_col is not in annotation dataframe columns
                # TypeError: annotation is None
                labels = None

            # This groupby iterates by row but we still need to grab the first
            # item because pandas doesn't cast the row to a vector
            pdf = os.path.join(tissue_method_path, row.index[0])

            i = None
            n = None
            if extra and '-of-' in extra:
                split = extra.split('-of-')
                i = int(split[0])
                n = int(split[-1])
                extra = None
            elif extra != '1-of-1':
                extra = None
            tex_generator = TeXGenerator(pdf, plottype, tissue,
                                        method, subset, groupby,
                                        i=i, n=n,
                                        labels=labels, extra=extra)

            # Weird hack for adding sections
            if prev_groupby != groupby or prev_subset != subset:
                j = 0
            else:
                j += 1

            print(f'\tsubset: {subset}, groupby: {groupby}, plottype: '
                  f'{plottype}, k: {k}, j: {j}, i: {i}, n: {n}, extra: {extra}')

            if j == 0 and tissue != 'Microbiome':
                tex += tex_generator.subsection_tex

                if subset.lower().startswith('subset'):
                    letter = subset.lower().split('subset')[-1].upper()
                    subset_col = 'subset' + letter
                    try:
                        subset_annotation = annotation.loc[annotation[subset_col]]
                    except KeyError:
                        sys.stderr.write(f'{subset_col} not found in {tissue} {method}!!!\n')
                        subset_annotation = annotation
                elif subset != 'allcells':
                    try:
                        subset_yaml = yaml_data['SUBSET'][groupby.upper()]
                        print(subset_yaml)
                        # Add table of counts
                        column = subset_yaml['FILTER_COLUMN']
                        value = subset_yaml['FILTER_VALUE']
                    except KeyError:
                        column = 'cell_ontology_class'
                        value = subset.lower().replace('_', ' ').rstrip('s')
                        if value.endswith('cell'):
                            value = value.split('cell')[0] + ' cell'
                    subset_annotation = annotation.query(f'{column} == "{value}"')
                else:
                    subset_annotation = annotation

                if 'expression' not in groupby_col:
                    try:
                        counts = subset_annotation.fillna("NA").groupby(groupby_col).size()
                        tex += tex_generator.make_table_tex(counts)
                    except KeyError:
                        # This is a custom plot and doesn't have a groupby
                        pass
                if groupby != 'highlighted':
                    # Don't count anything for highlighting cells
                    tex += r'''
\clearpage'''


#             else:
#                 tex += r'''
# \clearpage'''

            # Add section title and figure tex
            tex += tex_generator.figure_tex

            prev_subset = subset
            prev_groupby = groupby
            figures_added.append(pdf)

        # figures_not_added = parameters.index.difference(figures_added)
        # figures_not_added_str = '\n\t'.join(figures_not_added)
        # print(f'Figures not added: {figures_not_added_str}')

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

#!/usr/bin/env python3

# coding: utf-8

# This reads in a parameters file and generates an Rmd notebook
# for that organ and method.

# Usage: generate_from_template.py parameter_file.yaml
import sys
# import locale
#
# locale.setlocale(locale.LC_ALL, 'en_US.utf8')
# locale.setlocale(locale.LC_CTYPE, 'en_US.utf8')
import os

import yaml
import click

TISSUE = "TISSUE"
METHOD = "METHOD"
ADDITIONAL_CODE = 'ADDITIONAL_CODE'
SUBSET = "SUBSET"
FILTER_COLUMN = 'FILTER_COLUMN'
FILTER_VALUE = 'FILTER_VALUE'
RES = 'RES'
NPCS = "NPCS"
GENES = "GENES"
GROUPBY = "GROUPBY"
CODE_FOLDER = '29_tissue-specific_supplement_code'


def add_subset(name, filter_column, filter_value, res, npcs, genes, groupby):
    filter = f'rownames(tiss@meta.data)[grep("{filter_value}",tiss@meta.data${filter_column})]'
    genes_str = ', '.join(map(lambda x: f'"{x}"', genes))
    code = f"""{name}.cells.use = {filter}
{name}.n.pcs = {npcs}
{name}.res = {res}
{name}.genes_to_check = c({genes_str})
{name}.tiss = SubsetData(tiss, cells.use={name}.cells.use, )
{name}.tiss <- {name}.tiss %>% ScaleData() %>% 
  FindVariableGenes(do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5) %>%
  RunPCA(do.print = FALSE)
{name}.tiss <- {name}.tiss %>% FindClusters(reduction.type = "pca", dims.use = 1:sub.n.pcs, 
    resolution = sub.res.use, print.output = 0, save.SNN = TRUE) %>%
    RunTSNE(dims.use = 1:sub.n.pcs, seed.use = 10, perplexity=30)
"""
    if groupby is not None:
        # Append this subset's groupby to the list
        code += '''\n# Append this subset's groupby to the list
group.bys = c(group.bys, "{groupby}")'''

    code += '''dot_tsne_violin({name}.tiss, {name}.genes_to_check,
    save_folder, prefix = {name}, group.bys)
'''

    codeblock = f'''```{{r}}
{code}
```'''
    return codeblock




@click.command()
@click.argument('parameters_yaml')
@click.option('--template-file', default='Template.Rmd')
@click.option('--suffix', default='_template.Rmd')
def main(parameters_yaml, template_file='Template.Rmd',
         suffix='_template.Rmd'):
    # print command line arguments

    with open(parameters_yaml) as f:
        parameters = yaml.load(f)

    with open(template_file) as f:
        template = f.read()
        for parameter, value in parameters.items():
            if parameter == ADDITIONAL_CODE:
                # If ADDITIONAL_CODE is set to True
                if value:
                    basename = parameters[TISSUE] + "_" + parameters[METHOD] + '.Rmd'
                    filename = os.path.join(
                        '..', CODE_FOLDER, basename)
                    with open(filename) as g:
                        template += g.read()

            elif parameter == SUBSET:
                for name, kv in value.items():
                    # Make all the keys have lowercase names
                    kv = {k.lower(): v for k, v in kv.items()}
                    subset_code = add_subset(name, **kv)
                    import pdb; pdb.set_trace()
                    template += f'\n## Subset: {name}\n\n{subset_code}'
            elif value is not None:
                template = template.replace("{" + parameter + "}", str(value))
            else:
                template = template.replace("{" + parameter + "}", '')

        if ADDITIONAL_CODE not in parameters:
            template = template.replace('{' + ADDITIONAL_CODE + '}',
                                        '# No additional code')

    outfile = parameters[TISSUE] + "_" + parameters[METHOD] + suffix
    with open(outfile, 'w') as f:
        f.write(template)

if __name__ == "__main__":
    main()

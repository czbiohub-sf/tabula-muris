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

import yaml
import click

ADDITIONAL_CODE = 'ADDITIONAL_CODE'

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
                if value:

            elif value is not None:
                template = template.replace("{" + parameter + "}", str(value))
            else:
                template = template.replace("{" + parameter + "}", '')

        if ADDITIONAL_CODE not in parameters:
            template = template.replace('{' + ADDITIONAL_CODE + '}',
                                        '# No additional code')

    outfile = parameters['TISSUE'] + "_" + parameters['METHOD'] + suffix
    with open(outfile, 'w') as f:
        f.write(template)

if __name__ == "__main__":
    main()

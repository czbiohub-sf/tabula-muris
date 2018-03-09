# coding: utf-8

# This reads in a parameters file and generates an Rmd notebook
# for that organ and method.

# Usage: generate_from_template.py parameter_file.yaml
import sys
import yaml
import click

@click.command()
@click.argument('parameters_yaml')
@click.option('--template-file', default='Template.Rmd')
@click.option('--suffix', default='_template.Rmd')
def main(parameter_yaml, template_file='Template.Rmd',
         suffix='_template.Rmd'):
    # print command line arguments

    with open(parameter_yaml) as f:
        parameters = yaml.load(f)

    with open(template_file) as f:
        template = f.read()
        for parameter, value in parameters.items():
            if value is not None:
                template = template.replace("{" + parameter + "}", str(value))
            else:
                template = template.replace("{" + parameter + "}", '')

    outfile = parameters['TISSUE'] + "_" + parameters['METHOD'] + suffix
    with open(outfile, 'w') as f:
        f.write(template)

if __name__ == "__main__":
    main()

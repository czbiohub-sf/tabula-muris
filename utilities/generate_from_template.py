# coding: utf-8

# This reads in a parameters file and generates an Rmd notebook
# for that organ and method.

# Usage: generate_from_template.py parameter_file.yaml
import sys
import yaml

def main():
    # print command line arguments
    parameter_file = sys.argv[1]
    template_file = "Template.Rmd"

    with open(parameter_file) as f:
        parameters = yaml.load(f)

    with open(template_file) as f:
        template = f.read()
        for parameter, value in parameters.items():
            if value is not None:
                template = template.replace("{" + parameter + "}", str(value))
            else:
                template = template.replace("{" + parameter + "}", '')

    outfile = parameters['TISSUE'] + "_" + parameters['METHOD'] + "_template" + ".Rmd"
    with open(outfile, 'w') as f:
        f.write(template)

if __name__ == "__main__":
    main()

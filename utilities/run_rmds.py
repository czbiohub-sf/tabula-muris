import glob
import os
import subprocess

import click

"""
May need to set "locales" for Unicode because ASCII is stupid.

e.g. for a ENglish, US machine:

export LC_ALL=en_US.UTF-8 
export LC_ALL=en_US.UTF-8
"""

@click.command()
@click.option('--folder', default='.')
def main(folder):
    globber = os.path.join(folder, '*.Rmd')
    rmds = glob.iglob(globber)

    for rmd in rmds:
        subprocess.call(f"""echo \"rmarkdown::render('{rmd}', clean=TRUE)\" 
    | R --slave > {rmd}.out 2>{rmd}.err""")


if __name__ == "__main__":
    main()
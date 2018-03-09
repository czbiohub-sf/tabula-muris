#!/usr/bin/env python3

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
        click.echo(f'Starting {rmd} ...')
        command = ['echo', f'"rmarkdown::render(\'{rmd}\', clean=TRUE)"']
        # command = f"echo \"rmarkdown::render(\'{rmd}\', clean=TRUE)\" | R --slave"
        echo = subprocess.run(command, stdout=subprocess.PIPE)
        # echo = subprocess.Popen(command, stdout=subprocess.PIPE)
        rmd_output = subprocess.run(['R', '--slave'], stdin=echo.stdout,
                                    stdout=subprocess.PIPE)
        # echo.communicate()
        stdout = rmd + '.out'
        stderr = rmd + '.err'

        with open(stdout, 'w') as f:
            f.write(rmd_output.stdout)

        with open(stderr, 'w') as f:
            f.write(rmd_output.stderr)



    #     subprocess.call(f"""echo \"rmarkdown::render(\'{rmd}\', clean=TRUE)\"
    # | R --slave > {rmd}.out 2>{rmd}.err""")


if __name__ == "__main__":
    main()

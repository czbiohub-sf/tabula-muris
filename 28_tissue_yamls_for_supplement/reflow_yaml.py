import click
import yaml


def listify(genes):
    return genes.replace("'", '').replace(' ', '').split(',')

@click.command()
@click.argument('yamls', nargs=-1)
def cli(yamls):
    for filename in yamls:

        with open(filename) as f:
            data = yaml.load(f)

        data['GENES'] = listify(data["GENES"])
        if 'SUBSET' in data:
            for name, subset in data['SUBSET'].items():
                subset['GENES'] = listify(subset['GENES'])

        reflowed_data = yaml.dump(data,
                                  # Add the --- at the beginning of the file
                                  explicit_start=True,
                                  # Turn flow style to false to get
                                  # "bullet point" arrays
                                  default_flow_style=False)
        with open(filename + '.reflowed', 'w') as g:
            g.write(reflowed_data)


if __name__ == "__main__":
    cli()

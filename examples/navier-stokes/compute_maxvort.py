import click
import numpy as np

@click.command()
@click.argument('outfiles', nargs=-1, type=click.Path(exists=True))
def compute_maxvort(outfiles):

    for n, outfile in enumerate(sorted(outfiles)):
        sol = np.fromfile(outfile, dtype=np.float64)
        print n, abs(sol).max()

if __name__ == '__main__':
    compute_maxvort()

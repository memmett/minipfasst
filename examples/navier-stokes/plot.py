"""Parse output (from stdout) of NS runs and plot enstrophy."""

import click
import pylab
import numpy as np

from collections import namedtuple

pftuple = namedtuple('pftuple', ['step', 'iteration', 'level', 'data'])

@click.group()
def plot():
    pass

def parse_enstrophy(outfile):

    with open(outfile, 'r') as f:
        out = f.read()

    enstrophy = []
    for l in out.split('\n'):
        if 'enstrophy' in l:
            r = l.split()
            enstrophy.append((int(r[1]), float(r[2])))

    return sorted(enstrophy)


@plot.command()
@click.argument('outfiles', nargs=-1, type=click.Path(exists=True))
@click.option('--output', '-o', type=str, default=None)
@click.option('--start', '-s', type=int, multiple=True, default=[])
def enstrophy(outfiles, output, start):

    start = list(start) + len(outfiles) * [ 0 ]

    for i, outfile in enumerate(outfiles):
        enstrophy = parse_enstrophy(outfile)
        x, y = map(np.asarray, zip(*enstrophy))
        x = start[i] + x
        pylab.plot(x, y, label=outfile)

    pylab.legend(loc='upper left')
    pylab.ylim([0, 10.0])

    if output:
        pylab.savefig(output)
    else:
        pylab.show()


def parse_step_iter_level(outfile, variable):

    with open(outfile, 'r') as f:
        out = f.read()

    data = []
    for l in out.split('\n'):
        if variable in l:
            r = l.split()
            data.append(pftuple(int(r[1]), int(r[2]), int(r[3]), float(r[4])))

    return sorted(data)


pens = {
     (4, 1): { 'color': 'black', 'linestyle': '-', 'marker': '$1$', },
     (4, 2): { 'color': 'black', 'linestyle': '-', 'marker': '$2$', },
     (4, 3): { 'color': 'black', 'linestyle': '-', 'marker': '$3$', },
     (4, 4): { 'color': 'black', 'linestyle': '-', 'marker': '$4$', },
     (8, 1): { 'color': 'blue',  'linestyle': '-', 'marker': '$1$', },
     (8, 2): { 'color': 'blue',  'linestyle': '-', 'marker': '$2$', },
     (8, 3): { 'color': 'blue',  'linestyle': '-', 'marker': '$3$', },
     (8, 4): { 'color': 'blue',  'linestyle': '-', 'marker': '$4$', },
    (16, 1): { 'color': 'red',   'linestyle': '-', 'marker': '$1$', },
    (16, 2): { 'color': 'red',   'linestyle': '-', 'marker': '$2$', },
    (16, 3): { 'color': 'red',   'linestyle': '-', 'marker': '$3$', },
    (16, 4): { 'color': 'red',   'linestyle': '-', 'marker': '$4$', },
}

def nprocs(outfile):
    return int(outfile[4:6])


def plot_step_iter_level(outfiles, variable, logscale, output):

    for outfile in outfiles:
        pftuples = parse_step_iter_level(outfile, variable)
        maxlevel = max([ x.level for x in pftuples ])
        maxiteration = max([ x.iteration for x in pftuples ])
        for k in range(1, maxiteration + 1):
            xy = [ (x.step, x.data) for x in pftuples if x.level == maxlevel and x.iteration == k ]
            x, y = zip(*xy)
            if logscale:
                pylab.semilogy(x, y, label=outfile + ' iter ' + str(k), **pens[nprocs(outfile), k])
            else:
                pylab.plot(x, y, label=outfile + ' iter ' + str(k), **pens[nprocs(outfile), k])

    pylab.legend(loc='upper left')
    # pylab.ylim([0, 10.0])

    if output:
        pylab.savefig(output)
    else:
        pylab.show()


@plot.command()
@click.argument('outfiles', nargs=-1, type=click.Path(exists=True))
@click.option('--output', '-o', type=str, default=None)
def vorticity(outfiles, output):
    plot_step_iter_level(outfiles, 'maxvort', False, output)


@plot.command()
@click.argument('outfiles', nargs=-1, type=click.Path(exists=True))
@click.option('--output', '-o', type=str, default=None)
def residual(outfiles, output):
    plot_step_iter_level(outfiles, 'residual', True, output)


if __name__ == '__main__':
    plot()

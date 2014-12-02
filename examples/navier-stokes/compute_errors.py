"""Parse output (from stdout) of NS runs and plot enstrophy."""

import click
import numpy as np
import re
import os
import multiprocessing

from collections import namedtuple

pftuple = namedtuple('pftuple', ['step', 'iteration', 'level', 'data'])


def read(dname):
  """Read output directory *dname* and return list of available
  solutions.
  """

  prog = re.compile('(s(\d+)i(\d+)l(\d+)).dat')

  solutions = []
  for fname in os.listdir(dname):
    m = prog.search(fname)
    if m:
      step, iteration, level = map(int, m.groups()[1:])
      solutions.append(pftuple(step, iteration, level,
                               os.path.join(dname, m.group(1))))

  return solutions


def load(fname):
    return np.fromfile(fname + '.dat', dtype=np.float64)


def compute_error(t):
    s, i, ref, app = t
    xref = load(ref)
    xapp = load(app)
    return s, i, max(abs(xref-xapp))


@click.command()
@click.argument('reference', type=click.Path(exists=True))
@click.argument('approximate', type=click.Path(exists=True))
@click.option('--nprocs', '-n', default=1)
def compute_errors(reference, approximate, nprocs):

    ref = read(reference + '/out.d')
    app = read(approximate + '/out.d')

    maxiter = max([ x.iteration for x in ref ])
    refmap = { x.step: x.data for x in ref if x.iteration == maxiter }

    maxlev = max([ x.level for x in app ])
    appmap = { (x.step, x.iteration): x.data for x in app if x.level == maxlev }

    steps = set([ x.step for x in app ])
    iters = set([ x.iteration for x in app ])

    t = []
    for s in sorted(steps):
        for i in sorted(iters):
            try:
                t.append((s, i, refmap[s], appmap[s,i]))
            except:
                print 'WARNING: something funny with', s, i

    pool = multiprocessing.Pool(processes=nprocs)
    errors = pool.map(compute_error, t)
    pool.close()
    pool.join()

    for s, i, err in errors:
        print s, i, err


if __name__ == '__main__':
    compute_errors()

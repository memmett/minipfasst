"""Parse output (from stdout) of NS runs and plot enstrophy."""

import click
import re
import os
import multiprocessing
import subprocess

from collections import namedtuple


pftuple = namedtuple('pftuple', ['step', 'iteration', 'level', 'data'])


def read(dname):
  """Read output directory *dname* and return list of available
  solutions.
  """

  prog = re.compile('s(\d+)i(\d+)l(\d+).dat')

  solutions = []
  for fname in os.listdir(dname):
    m = prog.search(fname)
    if m:
      step, iteration, level = map(int, m.groups())
      solutions.append(pftuple(step, iteration, level,
                               os.path.join(dname, m.group(0))))

  return solutions


def compute_error(work):
  """Compute error"""
  ref, sols = work
  dname = os.path.dirname(os.path.abspath(__file__))
  out   = subprocess.check_output([os.path.join(dname, 'errors'), ref] + sols)

  errs = []
  prog = re.compile('s(\d+)i(\d+)l(\d+).dat\s*(\S+)')
  for l in out.split('\n'):
    m = prog.search(l)
    if m:
      step, iteration, level = map(int, m.groups()[0:3])
      error = float(m.group(4))
      errs.append((step, iteration, error))

  return errs


@click.command()
@click.argument('reference', type=click.Path(exists=True))
@click.argument('approximate', type=click.Path(exists=True))
@click.option('--nprocs', '-n', default=1)
def compute_errors(reference, approximate, nprocs):
  """Compute errors for each step and iteration."""

  ref = read(reference + '/out.d')
  app = read(approximate + '/out.d')

  maxiter = max([ x.iteration for x in ref ])
  refmap = { x.step: x.data for x in ref if x.iteration == maxiter }

  maxlev = max([ x.level for x in app ])
  appmap = { (x.step, x.iteration): x.data for x in app if x.level == maxlev }

  steps = sorted(refmap.keys())
  iters = sorted(set([ x.iteration for x in app ]))

  # build list of tuples: (step s ref. sol., [ step s sols. ])
  work = [ (refmap[s],
            [ appmap[s,i] for i in iters if (s,i) in appmap ]) for s in steps ]

  pool = multiprocessing.Pool(processes=nprocs)
  errors = pool.map(compute_error, work)
  pool.close()
  pool.join()

  for t in errors:
    for s, i, err in t:
      print s, i, err


if __name__ == '__main__':
  compute_errors()

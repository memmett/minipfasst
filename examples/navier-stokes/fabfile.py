"""Fabric (fabfile.org) tasks for PFASST paper."""

import glob
import io
import os
import subprocess

import numpy as np
import pandas

from fabric.api import *
from fabric.colors import *
from fabric.utils import *
from fabric.contrib.files import exists
from fabric.operations import get

from ggplot import *

env.minipfasst = '~/projects/minipfasst/'
env.scratch    = '/global/scratch2/sd/memmett/PFASST/tg/'

env.nprocs = [ 4, 8, 16 ]
env.niters = { 4: [4, 6],
               8: [4, 6, 8, 10],
               16: [4, 8, 12, 16, 20] }

def stage(trial, **kwargs):

    puts(green("staging trial " + trial))

    dname = os.path.join('stage.d', trial)
    local("mkdir -p " + dname)

    params = {
        # probin
        'nu': 0.01,
        'npts': 64,
        'niters': 4,
        'dt':  0.01,
        'nnodes': 3,
        'qtype': 1,
        'nthreads': 1,
        'nsteps': 8,
        'midpoint': 'f',
        'nlevs': 1,
        'input': '',

        # pbs
        'trial': trial,
        'queue': 'debug',
        'width': 1,
        'walltime': '00:30:00',
        }

    params.update(**kwargs)
    for param in [ 'nnodes', 'npts' ]:
        if isinstance(params[param], list):
            params[param] = ', '.join(map(str, params[param]))

    with open('probin.nml', 'r') as f:
        probin = f.read()
    probin = probin.format(**params)
    with open(os.path.join(dname, 'probin.nml'), 'w') as f:
        f.write(probin)

    with open('submit.qsub', 'r') as f:
        submit = f.read()
    params['mppwidth'] = 12 * params['width']
    submit = submit.format(**params)
    with open(os.path.join(dname, 'submit.qsub'), 'w') as f:
        f.write(submit)


def rsync():
    local('rsync -auz stage.d/ edison:{scratch}'.format(scratch=env.scratch))

@task
@hosts('edison.nersc.gov')
def build(clean=False):
    with cd(env.minipfasst + 'examples/navier-stokes'):
        run('git pull')
        if clean:
            run('make clean')
        run('make')

def trial_name(trial, nprocs=None, niters=None, nx=None, nlevs=2, reference=False, qtype='u', nnodes=3):
    if reference:
        tmpl = 's{trial}_ref_nx{nx:03d}'
    else:
        tmpl = 's{trial}_p{nprocs:02d}_k{niters:02d}_l{nlevs}_nx{nx:03d}_{qtype}{nnodes:d}'
    return tmpl.format(**locals())


@task
def taylor_green(trial, nx=128):
    # note that qtype=1028 is UNIFORM+NO_LEFT

    submit = []
    for p in env.nprocs:
        for k in env.niters[p]:
            name = trial_name(trial, nprocs=p, niters=k, nx=nx, qtype='gl', nnodes=5)
            stage(name, width=p,
                  dt=0.005, nsteps=64, nu='0.001d0', nlevs=2, npts=[64, 128], nnodes=[3, 5], qtype=1,
#                  nthreads=12, queue='regular', walltime='02:00:00', niters=k,
                  nthreads=12, queue='debug', walltime='00:30:00', niters=k,
                  input=env.scratch+'initial%snx%d.dat' % (trial, nx))
            submit.append("qsub {name}/submit.qsub".format(name=name))

    with open("stage.d/submit-all.sh", 'w') as f:
        f.write("#!/bin/sh\n")
        f.write("\n".join(submit))
        f.write("\n")

    rsync()


@task
def taylor_green_reference(trial, nx=256):

    name = trial_name(trial, nx=nx, reference=True)
    stage(name, width=1,
          dt=0.005, nsteps=64, nu='0.001d0', nlevs=1, npts=[nx], nnodes=[5], qtype=1,
          nthreads=12, queue='regular', walltime='16:00:00', niters=8,
          input=env.scratch+'initial%snx%d.dat' % (trial, nx))

    rsync()


@task
def compute_errors(trial):

    submit = """\
#!/bin/sh
#PBS -N comperrs
#PBS -q regular
#PBS -l mppwidth={width}
#PBS -l walltime={walltime}
#PBS -o /global/scratch2/sd/memmett/PFASST/tg/{trial}/errors
#PBS -e /global/scratch2/sd/memmett/PFASST/tg/{trial}/errors.stderr
#PBS -V

module load python
cd /global/scratch2/sd/memmett/PFASST/tg
python {minipfasst}/examples/navier-stokes/compute_errors.py -n {width} {reference} {trial}
"""

    with open('stage.d/submit-errors.sh', 'w') as all:
        all.write("#!/bin/sh\n")
        for p in env.nprocs:
            for k in env.niters[p]:
                tname = trial_name(trial, nprocs=p, niters=k, nx=128, qtype='gl', nnodes=5)
                rname = trial_name(trial, nx=256, reference=True)
                sname = '{trial}/submit-errors.qsub'.format(trial=tname)

                local('mkdir -p stage.d/' + tname)

                with open('stage.d/' + sname, 'w') as qsub:
                    qsub.write(submit.format(
                        trial=tname, width=6, walltime="00:20:00",
                        minipfasst=env.minipfasst, reference=rname))

                all.write("qsub {submit}\n".format(submit=sname))
    rsync()

@task
@hosts('edison.nersc.gov')
def plot_errors(trial, tname=None):

  if tname is None:
    puts(red("No trial name specified, please choose from one of:"))
    run("ls -d {scratch}/s{trial}*".format(scratch=env.scratch, trial=trial))
    return

  puts(yellow('generating error plot for: {trial}'.format(trial=tname)))

  get(remote_path='{scratch}/{trial}/errors'.format(scratch=env.scratch, trial=tname),
      local_path='{trial}.errors'.format(trial=tname))

  errors = subprocess.check_output("grep -r '^[0-9]' {trial}.errors".format(trial=tname), shell=True)
  df = pandas.read_table(io.BytesIO(errors), names=['step', 'iteration', 'error'], sep=' ')

  p = ggplot(df, aes(x='step', y='error', group='iteration', colour='factor(iteration)')) \
      + scale_y_log10() + geom_line()
  ggsave('{trial}.pdf'.format(trial=tname), p)
  puts(green('error plot saved in: {trial}.pdf'.format(trial=tname)))

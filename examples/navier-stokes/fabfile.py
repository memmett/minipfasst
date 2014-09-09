"""Fabric (fabfile.org) tasks for PFASST paper."""

import glob
import os.path

import numpy as np

from fabric.api import *
from fabric.colors import *
from fabric.utils import *
from fabric.contrib.files import exists

env.minipfasst = '~/projects/minipfasst/'
env.scratch    = '/global/scratch2/sd/memmett/PFASST/tg/'

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
        'nthreads': 1,
        'nsteps': 8,
        'midpoint': 'f',
        'nlevs': 1,

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

@task
def taylor_green():
    for p in [ 5, 10, 25, 50 ]:
        stage('p%02dl2nx128' % p, width=p,
              nsteps=1000, nu='0.001d0', nlevs=2, npts=[64,128], nnodes=[2,3],
              nthreads=8, queue='regular', walltime='02:00:00')
    rsync()

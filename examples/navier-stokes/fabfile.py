"""Fabric (fabfile.org) tasks for PFASST paper."""

import glob
import os.path

import numpy as np

from fabric.api import *
from fabric.colors import *
from fabric.utils import *
from fabric.contrib.files import exists


def NS(dname, exe='main.exe', probin='probin.nml', dry_run=False,
       **kwargs):

    puts(green("running trial " + dname))

    params = {
        'nu': 0.01,
        'npts': 64,
        'niters': 4,
        'dt':  0.01,
        'nnodes': 3,
        'nthreads': 1,
        'nsteps': 8,
        'midpoint': 'f',
        'nlevs': 1,
        }

    params.update(**kwargs)
    for param in [ 'nnodes', 'npts' ]:
        if isinstance(params[param], list):
            params[param] = ', '.join(map(str, params[param]))

    with open(probin, 'r') as f:
        probin = f.read()
    with open('/tmp/probin.nml', 'w') as f:
        f.write(probin.format(**params))

    run("mkdir -p " + os.path.join(dname, 'out.d'))
    run("cp " + exe + " " + dname)
    run("cp /tmp/probin.nml " + os.path.join(dname, 'probin.nml'))

    cmd = "cd {dname} && ./{exe} probin.nml | tee stdout 2> stderr".format(dname=dname,exe=os.path.basename(exe))

    if dry_run:
        puts(red(cmd))
    else:
        local(cmd)

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
        'walltime': '01:00:00',
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
    params['mppwidth'] = params['nthreads'] * params['width']
    submit = submit.format(**params)
    with open(os.path.join(dname, 'submit.qsub'), 'w') as f:
        f.write(submit)

def sync_stage():
    local('rsync -auz stage.d/ edison:/global/scratch2/sd/memmett/PFASST/tg/')


@task
def test(force_run=False):

    stage('p4l2nx128', nsteps=1000, nlevs=2, width=4, nthreads=8, npts=[64,128], nnodes=[2,3])

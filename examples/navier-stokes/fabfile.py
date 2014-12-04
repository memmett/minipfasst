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

env.nprocs = [ 4, 8, 16 ]
env.niters = { 4: [4, 6],
               8: [6, 8, 10],
               16: [8, 12, 16, 20] }

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

def trial_name(trial, nprocs, niters):
    return trial + '_p%02dk%02dl2nx256' % (nprocs, niters)

@task
def taylor_green(trial):
    # note that qtype=1028 is UNIFORM+NO_LEFT

    submit = []
    for p in env.nprocs:
        for k in env.niters[p]:
            name = trial_name(trial, p, k)
            stage(name, width=p,
                  dt=0.005, nsteps=64, nu='0.001d0', nlevs=2, npts=[64, 128], nnodes=[2, 3], qtype=1028,
                  nthreads=12, queue='regular', walltime='02:00:00', niters=k,
                  input=env.scratch+'initial%s.dat' % trial)
            submit.append("qsub {name}/submit.qsub".format(name=name))

    with open("stage.d/submit-all.sh", 'w') as f:
        f.write("#!/bin/sh\n")
        f.write("\n".join(submit))
        f.write("\n")

    rsync()


@task
def taylor_green_reference(trial):

    stage("reference%s" % trial, width=1,
          dt=0.005, nsteps=64, nu='0.001d0', nlevs=1, npts=[256], nnodes=[5], qtype=1028,
          nthreads=12, queue='regular', walltime='08:00:00', niters=8,
          input=env.scratch+'initial%s.dat' % trial)

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
                name = trial_name(trial, p, k)
                sname = '{trial}/submit-errors.qsub'.format(trial=name)

                with open('stage.d/' + sname, 'w') as qsub:
                    qsub.write(submit.format(
                        trial=name, width=6, walltime="01:00:00",
                        minipfasst=env.minipfasst, reference="reference"))

                all.write("qsub {submit}\n".format(submit=sname))
    rsync()


@task
def pull(trial):
    for p in [ 4, 8, 16 ]:
        name = trial + '_p%02dl2nx256' % p
        local("scp edison:{stdout} {lstdout}".format(
            stdout=env.scratch+name+"/stdout", lstdout=name + ".out"))

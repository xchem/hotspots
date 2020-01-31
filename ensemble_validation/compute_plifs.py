"""
Script to run  PLIP on all structures in an ensmeble.
This needs to be run in the PLIP environment (Python 2) on here:
/dls/science/groups/i04-1/software/anaconda/envs/PLIP

Write a bash script and execute it. Then, compile a csv dile that can be used as input for the PLIP grids
functions in ligand_grids.py
"""
from __future__ import print_function, division
from os.path import join, basename, dirname
import subprocess
from glob import glob
import shlex
import os
import stat

import luigi

def check_input(path):
    """
    PLIP doesn't take in .mol2 files, so if a .mol2 file is supplied, use openbabel to convert it to a pdb
    :param path: 
    :return: 
    """
    obabel_path = "/dls/science/groups/i04-1/software/anaconda/envs/PLIP/bin/obabel"
    extension = basename(path).split(".")[-1]

    if extension != 'pdb':
        print("converting mol2 to pdb")
        new_path = path.replace(extension, 'pdb')
        print(new_path)
        command = "{} {} -opdb -O {}".format(obabel_path, path, new_path)
        subprocess.call(shlex.split(command))

        return new_path

    else:
        return path

def write_plip_bash_script(ensemble_paths, ensemble_name):
    """
    
    :param ensemble_paths: 
    :return: 
    """
    fname = "run_{}_plips.sh".format(ensemble_name)
    plip_python_path = "/dls/science/groups/i04-1/software/anaconda/envs/PLIP/bin/python"

    out_paths = []
    with open(fname, "w") as f:
        f.write("#!/usr/bin/env bash \n")
        #f.write("alias plip='python /dls/science/groups/i04-1/software/mihaela/PLIP/plip/plip/plipcmd.py'")
        for p in ensemble_paths:
            in_path = check_input(p)
            out_path = dirname(in_path)
            f.write("{} /dls/science/groups/i04-1/software/mihaela/PLIP/plip/plip/plipcmd.py -f {} -xo {} &\n ".format(plip_python_path, in_path, out_path))
            out_paths.append(join(out_path, 'report.xml'))

    # Make the bash script executable by my user
    st = os.stat(fname)
    os.chmod(fname, st.st_mode | stat.S_IEXEC)

    command = "./{}".format(fname)
    subprocess.call(command)

    return out_paths

def compute_plip(ppath):
    plip_python_path = "/dls/science/groups/i04-1/software/anaconda/envs/PLIP/bin/python"

    in_path = check_input(path=ppath)
    out_path = dirname(in_path)
    command = "{} /dls/science/groups/i04-1/software/mihaela/PLIP/plip/plip/plipcmd.py -f {} -xo {}".format(plip_python_path, in_path, out_path)
    subprocess.call(shlex.split(command))

    return join(out_path, 'report.xml')

class lcompute_plip(luigi.Task):
    in_file = luigi.parameter.Parameter()

    def run(self):
        out_path =compute_plip(self.in_file)
        return out_path

    def output(self):
        out_file =join(dirname(self.in_file), 'report.xml')
        return luigi.LocalTarget(out_file)

if __name__ == "__main__":
    ensemble_paths =[x for x in  glob("/home/jin76872/Desktop/Mih/Data/ensemble_maps_validation/p38a/p38a/*/complex.mol2")]
    ensemble_name = "p38a"
    tar_f = write_plip_bash_script(ensemble_paths, ensemble_name)
    #p = check_input(ensemble_paths[0])



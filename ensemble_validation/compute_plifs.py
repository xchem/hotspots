from pathlib import Path
import subprocess
import shlex
import luigi

def check_input(path):
    """
    PLIP doesn't take non-pdb inputs, so if such a thing is  supplied, use openbabel to convert it to a pdb
    :param path: str or Path to input file.
    :return: 
    """
    # TODO: make into var
    obabel_path = "/dls/science/groups/i04-1/software/anaconda/envs/PLIP/bin/obabel"
    extension = Path(path).suffix

    if extension != 'pdb':
        new_path = path.replace(extension, 'pdb')
        print(new_path)
        command = "{} {} -opdb -O {}".format(obabel_path, path, new_path)
        subprocess.call(shlex.split(command))

        return new_path

    else:
        return path

def compute_plip(ppath):
    """
    
    :param ppath: 
    :return: 
    """
    # TODO: make into var. Also the PLIP executable in line 38.
    plip_python_path = "/dls/science/groups/i04-1/software/anaconda/envs/PLIP/bin/python"

    in_path = check_input(path=ppath).resolve()
    out_path = in_path.parent.resolve()
    command = "{} /dls/science/groups/i04-1/software/mihaela/PLIP/plip/plip/plipcmd.py -f {} -xo {}".format(plip_python_path, in_path, out_path)
    subprocess.call(shlex.split(command))

    return Path(out_path, 'report.xml')

class lcompute_plip(luigi.Task):
    in_file = luigi.parameter.Parameter()
    # TODO: add plip_python_path and obabel_path as parameters, which will be set up through the pipeline script

    def run(self):
        out_path =compute_plip(self.in_file)
        return out_path

    def output(self):
        out_file =Path(self.in_file.parent, 'report.xml')
        return luigi.LocalTarget(str(out_file.resolve()))

if __name__ == "__main__":
    ensemble_paths =[x for x in  Path("/home/jin76872/Desktop/Mih/Data/ensemble_maps_validation/p38a/p38a/").glob("*/complex.mol2")]
    ensemble_name = "p38a"
    #tar_f = write_plip_bash_script(ensemble_paths, ensemble_name)
    #p = check_input(ensemble_paths[0])



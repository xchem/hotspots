from __future__ import print_function, division
import json
from os.path import join
from Ensemble import EnsembleResult, MultipleEnsemble
from EnsembleReader import SIENAreader

def ensembleResultfromSIENA(ens_name, siena_dir, aligned=False):
    sr = SIENAreader(siena_dir=siena_dir,
                     ens_name= ens_name,
                     out_dir=join(root_dir, n))
    e = sr.output_ensemble()
    e.aligned = aligned
    return e


# Set some variables: TODO: import from json config file.

root_dir = "/home/jin76872/Desktop/Mih/Data/SIENA/"
ensemble_dict = {"BAZ2B": "4RVR",
                 "BRD1": "5N49"}
reference_ensemble = "BRD1"
SIENA_dirs = {key: join(root_dir, "{}_FILES".format(key)) for key in ensemble_dict.keys()}
aligned_references= False
hotspot_rotations = 3000
hotspot_charged_probes = False

########################### Run the pipeline ###########################

if len(ensemble_dict) == 1:
    n, s_dir = SIENA_dirs.items()[0]
    sr = SIENAreader(siena_dir=s_dir,
                     ens_name= ensemble_dict.keys,
                     out_dir=join(root_dir, n))
    e = sr.output_ensemble()
    e.aligned = aligned_references
    e.process_ensemble_proteins()
    e.make_ensemble_maps()

elif len(ensemble_dict) == 0:
    print("No ensembles given- type in the ensemble name in ensemble_dict")

else:
    e_dict = {key: ensembleResultfromSIENA(key, val) for key, val in SIENA_dirs.items()}
    me = MultipleEnsemble(ensemble_dict= e_dict,
                          ref_id=e_dict[reference_ensemble].ref_ID,
                          align_references=aligned_references)
    me.align_ensembles()
    me.binding_site_ensembles()


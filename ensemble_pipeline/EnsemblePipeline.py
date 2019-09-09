from __future__ import print_function, division
from SIENA_query import Search
import json
from os import mkdir
from os.path import join, exists
from Ensemble import EnsembleResult, MultipleEnsemble
from EnsembleReader import SIENAreader
from SIENA_query import Search
import tempfile
import pandas as pd

class EnsembleParams:
    def __init__(self, root_path):
        """
        json file containing the pipeline parameters
        :param str config_file: str, path to the pipeline's config file
        """
        self.pipeline_root = root_path
        self.config_file = join(root_path, "pipeline.json")
        # if config file has been supplied, load in the parameters from there
        if exists(self.config_file):
            with open(self.config_file, "r") as f:
                params =json.load(f)

            print("Existing config file found with parameters", params)
            self.ensemble_dict = params["ensemble_dict"]
            self.pipeline_root = params["pipeline_root"]
            self.reference_ensemble = params["reference_ensemble"]
            self.aligned_references = params["aligned_references"]
            self.hotspot_rotations = params["hotspot_rotations"]
            self.hotspot_charged_probes = params["hotspot_charged_probes"]
            self.SIENA_params = params["SIENA_params"]

        # if no config file present, create it with default values, tell the user what to do, and exit
        else:
            print("Could not find existsing file, v ELSA sum")
            param_dict = {"ensemble_dict": {"ensemble_name": ["PDB code", "ligand_name"]},
                          "pipeline_root": root_path,
                          "reference_ensemble": "",
                          "aligned_references": False,
                          "hotspot_rotations": 3000,
                          "hotspot_charged_probes": False,
                          "SIENA_params": {"reduction_procedure": "bb_clustering",
                                           "num_ensemble_members": str(50)}}
            p = json.dumps(param_dict, sort_keys=True, indent=4, separators=(',', ': '))

            with open(join(self.pipeline_root, "pipeline.json"), "w") as f:
                f.write(p)

class EnsembleIO:
    def __init__(self, root_path):
        files_file = join(root_path, "file_structure.json")
        if exists(files_file):
            with open(files_file, "r") as f:
                paths =json.load(f)

            self.SIENA_dirs = paths["SIENA_dirs"]
            self.ensemble_dirs = paths["ensemble_dirs"]
            self.ensemble_pdbs = paths["ensemble_pdbs"]
            self.multiple_ensembles_dirs = paths["multiple_ensembles_dirs"]
            self.hotspot_paths = paths["hotspot_paths"]
            self.binding_site_hotspot_paths = paths["binding_site_hotspot_paths"]
            self.ensemble_maps = paths["ensemble_maps"]

        else:
            self.SIENA_dirs = {}
            self.ensemble_dirs = {}
            self.ensemble_pdbs = {}
            self.multiple_ensembles_dirs = {}
            self.hotspot_paths = {}
            self.binding_site_hotspot_paths = {}
            self.ensemble_maps = {}

        self.params = EnsembleParams(root_path)
        self.ensembles = self.params.ensemble_dict.keys()

    def get_SIENA_dir(self, ensemble_name):
        # Check that the ensemble name is recognised:
        if ensemble_name not in self.ensembles:
            print("Unrecognised ensemble {}. Expected one of :{}".format(ensemble_name, self.ensembles))
            return
        # create the directory
        siena_dir = join(self.params.pipeline_root, "{}_SIENA".format(ensemble_name))
        if not exists(siena_dir):
            mkdir(siena_dir)
        else:
            print("Found existing SIENA directory at {}".format(siena_dir))

        self.SIENA_dirs[ensemble_name] = siena_dir
        return siena_dir

    def get_ensemble_dir(self, ensemble_name):
        if ensemble_name not in self.ensembles:
            print("Unrecognised ensemble {}. Expected one of :{}".format(ensemble_name, self.ensembles))
            return
        ens_dir = join(self.params.pipeline_root, ensemble_name)
        if not exists(ens_dir):
            mkdir(ens_dir)
        else:
            print("Found existing ensemble directory at {}".format(ens_dir))
        self.ensemble_dirs[ensemble_name] = ens_dir
        return ens_dir

    def get_hotspot_paths(self, ensemble_name):
        if len(self.ensemble_pdbs[ensemble_name]) == 0:
            print("No processed proteins found for ensemble {}. Update EnsembleIO or EnsembleResult.process()".format(ensemble_name))
            return
        else:
            hotspot_paths = [join(e, "fullsize_hotspots") for e in self.ensemble_pdbs[ensemble_name]]
            self.hotspot_paths = hotspot_paths
        return hotspot_paths

    def get_binding_site_hotspot_paths(self, ensemble_name):
        if len(self.ensemble_pdbs[ensemble_name]) == 0:
            print("No processed proteins found for ensemble {}. Update EnsembleIO or EnsembleResult.process()".format(ensemble_name))
            return
        else:
            binding_site_hotspot_paths = [join(e, "binding_site_hotspots") for e in self.ensemble_pdbs[ensemble_name]]
            self.binding_site_hotspot_paths = binding_site_hotspot_paths
        return  binding_site_hotspot_paths

    def update(self):
        param_dict = {key:val for key, val in self.__dict__.items() if key != "params"}
        p = json.dumps(param_dict, sort_keys=True, indent=4, separators=(',', ': '))
        with open(join(self.params.pipeline_root, "file_structure.json"), "w") as f:
            f.write(p)



class EnsemblePipeline:
    """
    Class that handles the assembly and alignment of ensembles
    """

    def __init__(self, pipeline_root):
        """
        json file containing the pipeline parameters
        :param str config_file: str, path to the pipeline's config file
        """
        self.params = EnsembleParams(pipeline_root)
        self.io = EnsembleIO(pipeline_root)

    def get_SIENA_ensemble(self, pdb, ligand):
        searcher = Search()
        ensemble = searcher.create_ensemble(pdb_code=pdb,
                                            ligand=ligand,
                                            reduction_procedure=self.params.SIENA_params["reduction_procedure"],
                                            num_members=str(self.params.SIENA_params["num_ensemble_members"]))

        return ensemble

    def ensembleResultfromSIENA(self, ens_name, siena_dir, ref_pdb):
        sr = SIENAreader(siena_dir=siena_dir,
                         ens_name=ens_name,
                         ref_pdb=ref_pdb,
                         out_dir= self.io.get_ensemble_dir(ens_name))
        e = sr.output_ensemble()
        e.aligned = self.params.aligned_references
        return e

    def run(self):
        # Run SIENA for the ensemble:
        SIENA_dict = {}
        for ens_name in self.params.ensemble_dict.keys():
            print(ens_name)
            pdb, lig = self.params.ensemble_dict[ens_name]
            print(type(str(pdb)), type(str(lig)))
            #siena_result = self.get_SIENA_ensemble(str(pdb), str(lig))
            s_dir = self.io.get_SIENA_dir(ens_name)
            #siena_result.save(out_dir=s_dir)
            print(s_dir)
            SIENA_dict[ens_name] = s_dir

        print(SIENA_dict)

        e_dict = {key: self.ensembleResultfromSIENA(key, val, self.params.ensemble_dict[key][0]) for key, val in SIENA_dict.items()}

        me = MultipleEnsemble(root_dir=self.params.pipeline_root,
                              ensemble_dict=e_dict,
                              ref_id=e_dict[self.params.reference_ensemble].reference_ID,
                              align_references=self.params.aligned_references)
        #me.align_references()
        e_pdbs = {}
        print("Aligning proteins within the ensembles")
        for e_name, ens in me.ensembles.items():
            #e_pdbs[e_name] = ens.process_ensemble_proteins()
            e_pdbs[e_name] = []
            e_dir = join(ens.root_dir)
            e_info = pd.read_csv(join(e_dir, "resultStatistic_fragments.csv"))
            prot_paths = []
            for idx, row in e_info.iterrows():
                p_name = "{}_{}-{}".format(e_name, row["PDB code"].strip(), row["PDB chains"].strip())
                p_path = join(e_dir, p_name, "{}.pdb".format(p_name))
                print(p_path)
                e_pdbs[e_name].append(p_path)
        self.io.ensemble_pdbs = e_pdbs
 
        print(me.ensembles.values())
        print("Calculating hotspots...")

        hot_path_dic = {e_name: ens.run_hotspots(paths_list= self.io.ensemble_pdbs[ens.ensemble_ID],
                                                          nrotations=int(self.params.hotspot_rotations),
                                                          charged=self.params.hotspot_charged_probes) for e_name, ens in me.ensembles.items()}

        self.io.hotspot_paths = hot_path_dic

        shrunk_hot_paths = {key: me.shrink_hotspots(hot_paths) for key, hot_paths in hot_path_dic.items()}
        self.io.binding_site_hotspot_paths = shrunk_hot_paths

        self.io.update()

if __name__ == "__main__":
    tempfile.tempdir = "/home/jin76872/Desktop/Mih/Data/tmp_superstar_ghecom"
    ep = EnsemblePipeline("/home/jin76872/Desktop/Mih/Data/ACS_fall_meeting_2019_slides")
    ep.run()




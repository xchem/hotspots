# Generic Python imports
from pathlib import Path
import logging
import tempfile

#Pipelining
import luigi


# Science Python imports
import pandas as pd

from run_hotspots_job import ParalleliselRunner
from compute_plifs import lcompute_plip

tempfile.tempdir =  "/home/jin76872/Desktop/Mih/Data/tmp_superstar_ghecom"

class Pipeline:
    def __init__(self, ensemble_name, ensemble_dir, aligned_protein_paths, hotspot_rotations=3000, hotspot_spheres=False, protonated=False, charged=False):
        """
        
        :param ensemble_name str: Name of the ensemble (eg 'p38a', 'p38a_out, etc') 
        :param ensemble_dir: str or Path to where pipeline output will be stored (apart from the PLIP result files, which are saved with the pdbs).
        :param aligned_protein_paths: str or Paths to where the aligned input pdbs live
        :param hotspot_rotations: probe rotations for the hotspots algorithm
        :param hotspot_spheres: whether to use the hotspot spheres or not
        :param protonated: bool, whether the input structures have been previously protonated (True) or if protonation is to be done through the ccdc API
        """
        self.name = ensemble_name
        self.ens_dir = ensemble_dir
        self.protein_paths = aligned_protein_paths
        self.hs_rotations = hotspot_rotations
        self.hs_spheres = hotspot_spheres
        self.hs_protonated = protonated
        self.hs_charged=charged
        self.out_file = Path(self.ens_dir, f'{self.name}.csv')

        log_msg = f" Running ensemble with parameters: \n " \
                  f"Name: {self.name} \n" \
                  f"Ensemble directory: {self.ens_dir} \n" \
                  f"Hotspot rotations: {self.hs_rotations} \n" \
                  f"Hotspot use_spheres: {self.hs_spheres} \n" \
                  f"Structures pre-protonated: {self.hs_protonated} \n" \
                  f"Ensemble output file: {str(self.out_file.resolve())}"
        # create a log file for the run
        logging.basicConfig(filename=Path(self.ens_dir, f'{self.name}.log'), filemode='w', level=logging.INFO)
        logging.info(log_msg)

    def check_path(self, ppath):
        if not Path(ppath).exists():
            logging.warning(f"Expected output file {ppath} does not exist")
            return

        if type(ppath) is str:
            return ppath
        else:
            return str(ppath.resolve())

    def run_pipeline(self):

        str_paths = [str(x.resolve()) for x in self.protein_paths]

        plip_paths = [lcompute_plip(p) for p in str_paths]

        luigi.build([ParalleliselRunner(str_paths, self.hs_rotations, self.hs_charged, self.hs_protonated, self.hs_spheres)],
                               local_scheduler=True,
                               workers=17)
        if self.hs_spheres:
            hs_paths = [Path(p.parent, "fullsize_hotspots_{}_spheres".format(self.hs_rotations), "out.zip") for p in self.protein_paths]
        else:
            hs_paths = [Path(p.parent, "fullsize_hotspots_{}".format(self.hs_rotations), "out.zip") for p in self.protein_paths]

        plip_paths = [self.check_path(Path(pp.parent, 'report.xml')) for pp in self.protein_paths]
        hs_paths = [self.check_path(hp) for hp in hs_paths]

        # Now put them into the output file.
        df = pd.DataFrame()
        df['protein_path'] = str_paths
        df['hotspot_result_path'] = hs_paths
        df['plip_path'] = plip_paths

        df.to_csv(self.out_file)

if __name__ == "__main__":
    from reduce_KLIFS_structures import get_protein_paths


    # ens_root_path = Path("/home/jin76872/Desktop/Mih/Data/ensemble_maps_validation/p38a/p38a/reduced_ensemble")
    # ens_paths = [x for x in ens_root_path.glob("*/*.pdb")]

    ens_root_path = Path("/home/jin76872/Desktop/Mih/Data/ensemble_maps_validation/CDK2/CDK2")
    root_df = pd.read_csv(Path(ens_root_path.parent, "removed_duplicates.csv"))
    ens_paths = get_protein_paths(root_df, ens_root_path, mode="complex")

    p = Pipeline(ensemble_name="CK2_no_duplicates",
                 ensemble_dir=ens_root_path,
                 aligned_protein_paths=ens_paths,
                 hotspot_rotations=3000,
                 hotspot_spheres=False,
                 protonated=True)

    p.run_pipeline()




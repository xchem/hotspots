from ccdc.protein import Protein
import os
from os.path import join, exists
from glob import glob
from hotspots.calculation import Runner
from hotspots import hs_io
from ccdc.io import MoleculeWriter
import shutil
import datetime


class DiamondRunner(object):
    """
    Class to handle individual hotspot mapping runs on Diamond.
    :param str stem_dir: Path to directory holding information for the target protein.
    :param str prot_name: Name of the target protein, eg "NUDT5".
    :param str event_id: Number of the PanDDA event calculated (eg "x114" in NUDT5-x114).
    :param list chains: Which chains of the protein to use for the hotspots calculation. Default = ["A"]..
    :param str keyword: Any additional identifier needed: eg "bound" or "ground" state of PanDDA event.
    :param str protein_path: default is none - script will find and prepare protein. If supplied, must be
                            prepared for hotspots calculation in advance. 
    """
    def __init__(self, stem_dir, prot_name, event_id, chains=None, keyword=None, protein_path=None, out_dir=None):
        self.stem = stem_dir
        self.protein_name = prot_name
        self.event_id = event_id
        if not chains:
            self.chains = ["A"]
        else:
            self.chains = chains
        self.keyword = keyword
        self.protein_path = protein_path
        self.out_dir = out_dir

    def prepare_protein(self):
        """
        Prepares the protein for the Hotspots Run.
        :param str prot_path: Path to .pdb file 
        :return 'ccdc.protein.Protein': protein that the calculation will be run on
        """
        prot = Protein.from_file(self.protein_path)
        for chain in prot.chains:
            print(chain.identifier)
            if chain.identifier not in self.chains:
                prot.remove_chain(chain.identifier)
        prot.remove_all_waters()
        for ligand in prot.ligands:
            prot.remove_ligand(ligand.identifier)
        prot.remove_all_metals()
        prot.add_hydrogens()
        return prot

    def extract_ligands(self):
        """
        Finds any interesting (non-solvent) bound molecules.
        For now: XChem structures have the ligand as LIG
        :param str prot_path: path to input pdb file 
        :return: 
        """
        if not self.protein_path:
            self.protein_path = self.find_protein()
        if not self.out_dir:
            self.out_dir = self.make_savedir()

        prot = Protein.from_file(self.protein_path)
        prot.detect_ligand_bonds()
        ligs = [l for l in prot.ligands if "LIG" in l.identifier]

        # Look for which chains the ligands come from:


        with MoleculeWriter(join(self.out_dir, "ligands.mol2")) as writer:
            for ligand in ligs:
                writer.write(ligand)

    def find_protein(self):
        """
        Tries to find the target .pdb file
        :return: a "ccdc.protein.Protein" instance
        """
        target_paths = glob(join(self.stem, "*.pdb"))
        if len(target_paths) == 0:
            target_paths = glob(join(self.stem, "*", "*.pdb"))

        if len(target_paths)== 1 and self.event_id in target_paths[0]:
            return target_paths[0]

        else:
            if self.keyword:
                prot_paths = [p for p in target_paths if (self.keyword in p and self.event_id in p)]
            else:
                prot_paths = [p for p in target_paths if self.event_id in p]
            assert len(prot_paths)==1, "Cannot find single input structure - check input directory"
            return prot_paths[0]

    def make_savedir(self):
        """
        Creates directory where hotspot result will be saved
        :return: str, path to saved directory
        """
        if not exists(join(self.stem, "hotspot_results")):
            os.mkdir(join(self.stem, "hotspot_results"))

        str_chains = "".join(self.chains)
        if self.keyword:
            out = join(self.stem, "hotspot_results",
                       "{}_{}_{}_{}_results".format(self.protein_name, self.event_id, str_chains, self.keyword))
        else:
            out = join(self.stem, "hotspot_results", "{}_{}_{}_results".format(self.protein_name, self.event_id, str_chains))

        if not exists(out):
            os.mkdir(out)
        return out

    def _save_superstar_grids(self, hs_runner):
        """
        Saves and Zips the SuperStar grids from the hotspots.calculation.Runner
        :param hs_runner: 
        :return: 
        """
        ss_dir = join(self.out_dir, "superstar_grids")
        if not exists(ss_dir):
            os.mkdir(ss_dir)
        for s in hs_runner.superstar_grids:
            s.grid.write(join(ss_dir, "superstar_{}.ccp4".format(s.identifier)))
        shutil.make_archive(ss_dir, 'zip', ss_dir)
        shutil.rmtree(ss_dir)

    def log_runner(self, nrot=None):
        """
        Saves the settings of the class in a log file in 
        :return: 
        """
        if not self.out_dir:
            self.out_dir = self.make_savedir()

        if not self.protein_path:
            self.protein_path = self.find_protein()

        now = datetime.datetime.now()

        with open(join(self.out_dir, "hotspot_results.log"), "a") as f:
            f.write("Hotspots log calculated at {}-{}-{} {}:{} \n".format(now.year, now.month, now.day, now.hour, now.minute))
            for key, val in self.__dict__.items():
                f.write(str(key) + ": " + str(val) + "\n")
            if nrot:
                f.write("Number of rotations: {} \n".format(str(nrot)))
            f.write("\n")


    def run_hotspot_calculation(self, nrot=100000, method="ghecom", charged=True, sphere_maps=False, save_ligand=True):
        """
        Runs the hotspots calculation on the specified PDB structure
        :return: 
        """

        if not self.out_dir:
            self.out_dir = self.make_savedir()

        if not self.protein_path:
            self.protein_path = self.find_protein()
            protein = self.prepare_protein()
        else:
            protein = Protein.from_file(self.protein_path)

        if save_ligand:
            self.extract_ligands()

        # log the run parameters
        self.log_runner(nrot)

        h = Runner()
        settings = h.Settings()
        settings.nrotations = nrot
        settings.apolar_translation_threshold = 15
        settings.polar_translation_threshold = 15
        settings.sphere_maps = sphere_maps

        result = h.from_protein(protein=protein,
                                charged_probes=charged,
                                probe_size=7,
                                buriedness_method=method,
                                cavities=None,
                                nprocesses=5,
                                settings=settings)
        #self.out_dir = self.make_savedir()
        # Save and zip the SuperStar Grids:
        self._save_superstar_grids(h)

        # Save and zip the Results
        with hs_io.HotspotWriter(self.out_dir, visualisation="pymol", grid_extension=".ccp4", zip_results=True) as writer:
            writer.write(result)












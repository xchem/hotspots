from __future__ import print_function, division
from ccdc.protein import Protein
from ccdc.io import MoleculeWriter
from os.path import join, basename, dirname, exists
from os import mkdir
from hotspots.calculation import Runner
from hotspots import hs_io
from shutil import make_archive, rmtree
import sys
import luigi

class RunHotspots:
    def __init__(self, pdb_file, number_rotations, charged=False):
        self.in_file = pdb_file
        self.number_rotations = number_rotations
        self.charged = charged
        self.out_dir = join(dirname(self.in_file), "fullsize_hotspots")
        if not exists(self.out_dir):
            mkdir(self.out_dir)

    def prepare_protein(self):

        prot = Protein.from_file(self.in_file)
        prot.remove_all_waters()
        for ligand in prot.ligands:
            prot.remove_ligand(ligand.identifier)
        prot.remove_all_metals()
        prot.add_hydrogens()
        return prot

    def _save_superstar_grids(self, hs_runner):
        """
        Saves and Zips the SuperStar grids from the hotspots.calculation.Runner
        :param hs_runner: 
        :return: 
        """
        ss_dir = join(self.out_dir, "superstar_grids")
        if not exists(ss_dir):
            mkdir(ss_dir)
        for s in hs_runner.superstar_grids:
            s.grid.write(join(ss_dir, "superstar_{}.ccp4".format(s.identifier)))
        make_archive(ss_dir, 'zip', ss_dir)
        rmtree(ss_dir)

    def run_hotspot_calculation(self,  method="ghecom", sphere_maps=False):
        """
        Runs the hotspots calculation on the specified PDB structure
        :return: 
        """
        h = Runner()
        settings = h.Settings()
        settings.nrotations = self.number_rotations
        settings.apolar_translation_threshold = 15
        settings.polar_translation_threshold = 15
        settings.sphere_maps = sphere_maps

        result = h.from_protein(protein=self.prepare_protein(),
                                charged_probes=self.charged,
                                probe_size=7,
                                buriedness_method=method,
                                cavities=None,
                                nprocesses=3,
                                settings=settings)
        #self.out_dir = self.make_savedir()
        # Save and zip the SuperStar Grids:
        self._save_superstar_grids(h)

        # Save and zip the Results
        with hs_io.HotspotWriter(self.out_dir, visualisation="pymol", grid_extension=".ccp4", zip_results=True) as writer:
            writer.write(result)



class lRunner(luigi.Task):
    in_pdb = luigi.parameter.Parameter()
    nrot = luigi.parameter.IntParameter()
    charged = luigi.parameter.BoolParameter()

    def run(self):
        RunHotspots(self.in_pdb, number_rotations=self.nrot, charged=self.charged).run_hotspot_calculation()

    def output(self):
        tar = join(dirname(self.in_pdb), "fullsize_hotspots", "out.zip")
        return luigi.LocalTarget(tar)


class ParalleliselRunner(luigi.WrapperTask):

    in_list = luigi.parameter.ListParameter()
    nrot = luigi.parameter.IntParameter()
    charged = luigi.parameter.BoolParameter()

    def requires(self):
        for k in self.in_list:
            yield lRunner(k, self.nrot, self.charged)


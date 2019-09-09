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
    def __init__(self, pdb_file, number_rotations, charged=False, data_source='SIENA'):
        self.in_file = pdb_file
        self.number_rotations = number_rotations
        self.charged = charged
        self.data_source = data_source # Can be 'SIENA', 'XChem' or 'KLIFS'
        self.out_dir = join(dirname(self.in_file), "fullsize_hotspots_{}".format(self.number_rotations))
        if not exists(self.out_dir):
            mkdir(self.out_dir)

    def prepare_protein(self):

        prot = Protein.from_file(self.in_file)

        if self.data_source == 'SIENA':
            chains = list(basename(self.in_file).split("-")[1].replace(".pdb", ""))
            print(chains, prot.chains)
            for ch in prot.chains:
                if ch.identifier not in chains:
                    prot.remove_chain(ch.identifier)

        prot.remove_all_waters()
        for ligand in prot.ligands:
            prot.remove_ligand(ligand.identifier)
        prot.remove_all_metals()
        if self.data_source != 'KLIFS':
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
                                nprocesses=1,
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
    data_source = luigi.parameter.Parameter()

    def run(self):
        RunHotspots(self.in_pdb, number_rotations=self.nrot, charged=self.charged, data_source=self.data_source).run_hotspot_calculation()
        #RunHotspots(self.in_pdb, number_rotations=self.nrot, charged=self.charged).run_hotspot_calculation()

    def output(self):
        tar = join(dirname(self.in_pdb), "fullsize_hotspots_{}".format(self.nrot), "out.zip")
        #tar = join(dirname(self.in_pdb), "fullsize_hotspots_{}".format(self.nrot), "out")
        return luigi.LocalTarget(tar)


class ParalleliselRunner(luigi.WrapperTask):

    in_list = luigi.parameter.ListParameter()
    nrot = luigi.parameter.IntParameter()
    charged = luigi.parameter.BoolParameter()
    data_source = luigi.parameter.Parameter()

    def requires(self):
        for k in self.in_list:
            yield lRunner(k, self.nrot, self.charged, self.data_source)
            #yield lRunner(k, self.nrot, self.charged)

if __name__ == "__main__":
    import pandas as pd
    import tempfile
    tempfile.tempdir = "/home/jin76872/Desktop/Mih/Data/tmp_superstar_ghecom"
    from glob import glob
    from shutil import copy

    tar_dir = "/home/jin76872/Desktop/Mih/Data/HAO-1"
    ensembles = ["HAO1"]
    nrotations = 3000
    charged = False

    for e in ensembles:
    #     e_dir = join(tar_dir, e)
    #     e_info = pd.read_csv(join(e_dir, "resultStatistic_fragments.csv"))
    #     prot_paths = []
    #     for idx, row in e_info.iterrows():
    #         p_name = "{}_{}-{}".format(e, row["PDB code"].strip(), row["PDB chains"].strip())
    #         p_path = join(e_dir, p_name, "{}.pdb".format(p_name))
    #         print(p_path)
    #         prot_paths.append(p_path)
        prot_paths = glob(join(tar_dir,"*",  "*.pdb"))
        # new_paths = [join(dirname(p), basename(p).replace(".pdb", "")) for p in prot_paths]
        #
        # new_prot_paths = []
        #
        # for new, old in zip(new_paths, prot_paths):
        #     if not exists(new):
        #         mkdir(new)
        #     newpath = copy(old, new)
        #     new_prot_paths.append(newpath)

        luigi.build([ParalleliselRunner(prot_paths, nrotations, charged, data_source='other')], local_scheduler=True, workers=30)


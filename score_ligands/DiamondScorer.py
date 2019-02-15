from ccdc.protein import Protein
import os
from os.path import join, exists
from glob import glob
from hotspots import hs_io
from ccdc.io import MoleculeWriter, MoleculeReader
import numpy as np
import shutil


class DiamondScorer(object):
    """
    """
    def __init__(self, stem, prot_name, event_id, keyword=None, out_dir=None, chains=None):
        self.stem = stem
        self.keyword = keyword
        self.data_dir = None
        self.protein_name = prot_name
        self.event_id = event_id
        self.out_dir = out_dir
        self.chains = chains
        if not self.chains:
            self.chains = ["A"]



    def get_data_dir(self):
        res_dir = join(self.stem, "hotspot_results")
        assert exists(res_dir), "Cannot find results directory for the target"
        str_chains = "".join(self.chains)
        if self.keyword:
            out = join(res_dir, "{}_{}_{}_{}_results".format(self.protein_name, self.event_id, str_chains, self.keyword))

        else:
            out = join(self.stem, "hotspot_results", "{}_{}_{}_results".format(self.protein_name, self.event_id, str_chains))
            if not exists(out):
                print("Can't find event-specific results dir")

        self.data_dir = out

    def get_hotspot(self):
        """
        
        :return: 
        """
        if not self.data_dir:
            self.get_data_dir()

        hs_reader = hs_io.HotspotReader(join(self.data_dir, "out.zip"))
        hs_result = hs_reader.read()
        return hs_result

    def get_ligands(self):
        """
        
        :return: 
        """
        if not self.data_dir:
            self.get_data_dir()
        lig_reader = MoleculeReader(join(self.data_dir, "ligands.mol2"))

        return lig_reader

    def read_protein(self):
        """
        Reads in the protein used to calculate the hotspot maps in data_dir.
        :return: 
        """
        hr = self.get_hotspot()

        return hr.protein


    def get_ligand_chain(self, ligand):
        """
        For each ligand, determine which protein chain (or chains) it is closest to.
        :param ccdc.molecule.Molecule ligand: 
        :return: 
        """
        prot = self.read_protein()
        bs = Protein.BindingSiteFromMolecule(prot, ligand, 5)
        ligand_chains = list(set([r.identifier.split(":")[0] for r in bs.residues]))
        return ligand_chains


    def score_ligands(self):
        """
        
        :return: 
        """
        hs = self.get_hotspot()
        ligs = self.get_ligands()

        if len(ligs) == 0:
            print("No ligands to score")
            return

        if not self.out_dir:
            self.out_dir = self.data_dir

        scored_ligs = []

        for lig in ligs:
            scored_lig = hs.score(lig)
            # change the identifier of the scored_lig to include the protein chain it came from
            lig_chains = "".join(self.get_ligand_chain(scored_lig))
            ligand_score = np.mean([a.partial_charge for a in scored_lig.heavy_atoms])
            scored_lig.identifier += "_{}_{}".format(lig_chains, round(ligand_score, 2))
            scored_ligs.append(scored_lig)

        with MoleculeWriter(join(self.out_dir, "scored_ligands.mol2")) as writer:
            for ligand in scored_ligs:
                writer.write(ligand)

        return scored_ligs



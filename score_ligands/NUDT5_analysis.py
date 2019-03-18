from __future__ import print_function, division
from ccdc.protein import Protein
from ccdc.molecule import Molecule
import os
from os.path import join, exists, dirname
from glob import glob
from hotspots import hs_io
from ccdc.io import MoleculeWriter, MoleculeReader
import numpy as np
import shutil
import pandas as pd
from DiamondScorer import DiamondScorer

"""
Classes to handle downstream processing from DiamondRunner and DiamondScorer
"""


class DFragScoreSummary(object):
    """
    Class to handle working with all the scored fragments for a target.
    """
    def __init__(self, results_dir, target_name, output_dir=None):
        self.stem = results_dir
        self.prot_name = target_name
        self.out_dir = output_dir
        self.summary_DF = None

    def make_output_dir(self):
        """
        Makes an output directory if not already set.
        :return: 
        """
        if not self.out_dir:
            out = join(self.stem, "scored_fragments_summary_{}".format(self.prot_name))
            if not exists(out):
                os.mkdir(out)
            self.out_dir = out
        else:
            print("Output directory already set")

    def _read_ligands(self):
        """
        Reads the scored ligands from the result directories for the target.
        :return: a :class: ccdc.io.MoleculeReader instance
        """
        lig_paths = glob(join(self.stem, "*", "scored_ligands.mol2"))
        scored_ligs = MoleculeReader(lig_paths)
        return scored_ligs

    def _make_lig_dict(self, mol):
        """
        Mole has to be a scored molecule output by DiamondScorer.
        :param mol: a :class: ccdc.molecule.Molecule instance
        :return: python dictionary
        """
        mol_fields = mol.identifier.split("_")
        print(mol_fields)
        mol_dic = {"Target": mol_fields[1],
                   "Identifier": mol.identifier,
                   "Crystal_ID": mol_fields[2],
                   "Chains": mol_fields[3],
                   "Hotspot_Score": float(mol_fields[4]),
                   "Perc_Unfulfilled_Polar": self.get_percent_unfulfilled_polar(mol),
                   "Unfulfilled_Polar": [a.label for a in self.get_unfulfilled_polar_contacts(mol)],
                   "HLE": self.get_hotspot_ligand_efficiency(mol),
                   "Perc_Unscored_Atoms" :self.get_perc_unscored_atoms(mol),
                   "Unscored_Atoms": [a for a in mol.heavy_atoms if a.partial_charge==0.0],
                   "Smiles": mol.smiles}
        return mol_dic

    def make_summary_DF(self):
        """
        Create Pandas dataframe to hold all data for the scored ligands
        :return: 
        """
        sc_ligs = self._read_ligands()
        df_cols = ["Target", "Identifier", "Crystal_ID", "Chains", "Hotspot_Score", "Perc_Unfulfilled_Polar", "Unfulfilled_Polar", "HLE", "Perc_Unscored_Atoms", "Unscored_Atoms","Smiles"]
        sum_DF = pd.DataFrame(columns=df_cols)

        for lig in sc_ligs:
            sum_DF.loc[len(sum_DF)] = self._make_lig_dict(lig)

        self.summary_DF = sum_DF

    def save_summary_DF(self, format="csv"):
        """
        
        :param format: 
        :return: 
        """
        if not self.summary_DF:
            self.make_summary_DF()

        if not self.out_dir:
            self.make_output_dir()

        if format=="csv":
            self.summary_DF.to_csv(join(self.out_dir, "scored_fragments_summary.csv"))
        else:
            print("unsupported format")

    def get_percentile_fragments(self, percentile=90, mode="above", save=True):
        """
        Selects the top scoring fragments
        :param percentile: 
        :return: 
        """
        scores = self.summary_DF["Hotspot_Score"]
        threshold = np.percentile(scores.values, percentile)

        if mode== "above":
            select = self.summary_DF.loc[scores > threshold]
        elif mode == "below":
            select = self.summary_DF.loc[scores < threshold]
        else:
            print("Unrecognised value for param 'mode' for method DFragScoreSummary.get_percentile_fragments()")
            return

        ids = select["Identifier"].values
        sc_ligs = self._read_ligands()
        top_ligs = [lig for lig in sc_ligs if lig.identifier in ids]

        if save:
            with MoleculeWriter(join(self.out_dir, "fragments_{}_{}_percentile.mol2".format(mode, percentile))) as writer:
                for l in top_ligs:
                    writer.write(l)
        return top_ligs

    def get_feat_dic(self, feat):
        """
        Gets the dictionary of feature attributes
        :param feat: a :class: DiamondScorer.LigandFeature instance
        :return: python dictionary
        """
        feat_dic = {"Coordinates": feat.centroids,
                    "Ligand_Atom": feat.atom.label,
                    "Score": feat.atom_score,
                    "Distance_to_closest_point":feat.distances,
                    "Near_Residues_to_Centroid": feat.centroid_proxy_residues,
                    "Near_Residues_to_Atom": feat.scored_atom_proxy_residues,
                    "Type": feat.feature_type,
                    "Size": feat.feature_size}
        return feat_dic

    def get_preferred_features(self, tar_ligs, keyw="bound"):
        """
        Extracts features for the ligands in the list. 
        Summarises them.
        :return: 
        """

        if len(tar_ligs) == 0:
            print("No ligands found")
            return

        new_DF_cols = ["Ligand_ID", "Type", "Ligand_Atom", "Size", "Score", "Distance_to_closest_point", "Coordinates", "Near_Residues_to_Centroid", "Near_Residues_to_Atom"]
        new_DF = pd.DataFrame(columns=new_DF_cols)

        for t in tar_ligs:
            mol_dic = self._make_lig_dict(t)
            dsc = DiamondScorer(stem=res_dir,
                                prot_name=mol_dic["Target"],
                                x_id=mol_dic["Crystal_ID"],
                                keyword=keyw,
                                chains=list(mol_dic["Chains"]))
            feats = dsc.get_ligand_features(t)

            for f in feats:
                f_dic = self.get_feat_dic(f)
                f_dic["Ligand_ID"] = mol_dic["Identifier"]
                new_DF.loc[len(new_DF)] = f_dic

        return new_DF

    def get_unfulfilled_polar_contacts(self, mol):
        """
        Looks for donor and acceptor atoms in the ligand with scores of zero
        :param mol: a :class: ccdc.molecule.Molecule instance with assigned Hotspot scores in partial charge field.
        :return: list
        """
        unf_polar = [a for a in mol.atoms if (a.is_donor or a.is_acceptor) and a.partial_charge == 0]
        return unf_polar

    def get_percent_unfulfilled_polar(self, mol):
        """
        
        :param mol: 
        :return: 
        """
        unf = self.get_unfulfilled_polar_contacts(mol)
        all_polar = [m for m in mol.atoms if (m.is_donor or m.is_acceptor)]
        perc = len(unf)/len(all_polar)*100

        return perc

    def get_hotspot_ligand_efficiency(self, mol):
        """
        HLE = hotspot_score/(number of nonzero scored carbons). 
        :param mol: 
        :return: 
        """
        try:
            hs_score = np.mean([a.partial_charge for a in mol.heavy_atoms])
        except TypeError:
            hs_score = 0
        scored_carbons = [m for m in mol.heavy_atoms if (m.partial_charge > 0.0 and m.atomic_symbol=="C")]
        #print(scored_carbons)
        if len(scored_carbons) > 0:
            HLE = hs_score/len(scored_carbons)
        else:
            HLE=None
        return HLE

    def get_perc_unscored_atoms(self, mol):
        """
        Returns a list of atoms that have not been scored.
        :param mol: 
        :return: 
        """
        unscored = [m for m in mol.heavy_atoms if m.partial_charge == 0.0]
        perc = len(unscored)/len(mol.heavy_atoms)*100

        return perc

    def create_scorer(self, lig_ID):
        """
        TODO: Move into DiamondScorer.Redundant with get_prefereed_features
        Takes a ligand ID and returns a breakdown of how atoms are scored.
        :param str lig_ID: the assigned ligand ID (found in self.summary_DF["Ligand_ID"])
        :return: 
        """
        ligs = self._read_ligands()
        liglist = [l for l in ligs if l.identifier == lig_ID]

        assert len(liglist) == 1

        l = liglist[0]
        mol_dic = self._make_lig_dict(l)
        dsc = DiamondScorer(stem=res_dir,
                            prot_name=mol_dic["Target"],
                            x_id=mol_dic["Crystal_ID"],
                            keyword="bound",
                            chains=list(mol_dic["Chains"]))
        return dsc

if __name__ == "__main__":

    res_dir = "/home/jin76872/Desktop/Mih/Data/NUDT5A/NUDT5A_fragment_hits/hotspot_results"
    prot = "NUDT5"
    fss = DFragScoreSummary(res_dir, prot)
    #fss.save_summary_DF()
    fss.make_summary_DF()
    """
    #top_90 = fss.get_percentile_fragments(percentile=10, mode="below", save=False)
    l = [n for n in fss._read_ligands() if n.identifier=="F:LIG1_NUDT5_0681_AB_19.33"][0]
    mol_dic = fss._make_lig_dict(l)
    dsc = DiamondScorer(stem=res_dir,
                        prot_name=mol_dic["Target"],
                        x_id=mol_dic["Crystal_ID"],
                        keyword="bound",
                        chains=list(mol_dic["Chains"]))
    features = dsc.get_ligand_features(l)

    feat_df = fss.get_preferred_features([l])
    """






    













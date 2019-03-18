from __future__ import print_function, division
from ccdc.protein import Protein
from ccdc.molecule import Molecule
import os
from os.path import join, exists, dirname
from glob import glob
from ccdc.io import MoleculeWriter, MoleculeReader, _CSDDatabaseLocator
from ccdc.docking import Docker
from DiamondScorer import DiamondScorer
import numpy as np

from ccdc.utilities import _private_importer
with _private_importer():
    import ChemicalAnalysisLib
    import ConformerGeneratorLib


class FollowUp(object):
    """
    Class to handle the docking and scoring of follow-ups, assuming the binding mode of the initial fragment is retained.
    """

    def __init__(self, followup_path, hotspot_path, fragment=None):
        """

        :param str followup_path: path to csv file with follow-up Smiles
        :param str hotspot_path: path to where the hotspot that is used for scoring will be.
        :param fragment: a :class: ccdc.molecule.Molecule instance
        """
        self.csv_path = followup_path
        self.hotspot_path = hotspot_path
        self.reference_fragment = fragment

    def from_smiles(self, smiles, identifier=None, generate_initial_sites=True):
        """
        Function taken from Pete's code.
        :param smiles: 
        :param identifier: 
        :param generate_initial_sites: 
        :return: 
        """
        if identifier is None:
            identifier = smiles

        if generate_initial_sites:
            parameter_files = _CSDDatabaseLocator.get_conformer_parameter_file_location()
            molmaker = ConformerGeneratorLib.MoleculeTo3D(parameter_files)
            mol = Molecule(identifier, molmaker.create_conformation(smiles))
        else:
            molmaker = ChemicalAnalysisLib.SMILESMoleculeMaker()
            mol = Molecule(identifier, _molecule=molmaker.siteless_atoms(smiles))
        return mol

    def read_followups(self):
        """
        Reads in the follow-ups from a .csv file
        :return: list of ccdc.molecule.Molecule instances
        """
        mols = []
        with open(self.csv_path, "r") as f:
            for line in f.readlines():
                print(line)
                if len(line)>1:
                    mol = self.from_smiles(line)
                    mols.append(mol)

        with MoleculeWriter(os.path.join(self.hotspot_path, "follow_ups.mol2")) as writer:
            for ligand in mols:
                writer.write(ligand)


    def get_fragment(self):
        """
        Gets the reference fragment
        :return: 
        """
        if not self.reference_fragment:
            ligs_p = join(self.hotspot_path, "scored_ligands.mol2")
            ligs = MoleculeReader(ligs_p)
            self.reference_fragment = ligs[0]

    def get_scorer_result(self):
        """
        
        :return: 
        """
        if not self.reference_fragment:
            self.get_fragment()

        mol_fields = self.reference_fragment.identifier.split("_")
        print(mol_fields)
        mol_dic = {"Target": mol_fields[1],
                   "Identifier": self.reference_fragment.identifier,
                   "Crystal_ID": mol_fields[2],
                   "Chains": mol_fields[3],
                   "Hotspot_Score": float(mol_fields[4])}

        dsc = DiamondScorer(stem= dirname(self.hotspot_path),
                            prot_name=mol_dic["Target"],
                            x_id=mol_dic["Crystal_ID"],
                            keyword="bound",
                            chains=list(mol_dic["Chains"]))
        return dsc


    def run_docking(self):
        """
        Reads in the follow-ups and tries to dock them
        :return: 
        """
        self.get_fragment()
        self.read_followups()

        docker = Docker()
        settings = docker.settings
        tempd = join(self.hotspot_path, "docking_tmp")

        # Get out the reference protein:
        scorer = self.get_scorer_result()
        hs = scorer.get_hotspot()
        prot = hs.protein

        # Change this from DiamondRunner - to save the protein in the results directory
        with MoleculeWriter(join(self.hotspot_path, "protein.pdb")) as prot_writer:
            prot_writer.write(prot)

        settings.add_protein_file(join(self.hotspot_path, "protein.pdb"))
        settings.binding_site = settings.BindingSiteFromPoint(settings.proteins[0], self.reference_fragment.centre_of_geometry(), 10.0)
        settings.fitness_function = 'plp'
        settings.autoscale = 10.0
        settings.output_directory = tempd
        #settings.output_directory = self.in_dir
        settings.output_file = "docked_ligands.mol2"
        settings.add_ligand_file(join(self.hotspot_path, "follow_ups.mol2"), ndocks = 10)

        # setup constraints
        settings.add_constraint(settings.TemplateSimilarityConstraint(type="all", template=self.reference_fragment, weight=150))
        results = docker.dock()
        output_file = os.path.join(settings.output_directory, settings.output_file)
        docked_molecules = [m for m in MoleculeReader(os.path.join(tempd, output_file))]

        return docked_molecules

    def score_docked_ligands(self):
        #docked_mols = self.run_docking()
        docked_mols = MoleculeReader(join(self.hotspot_path,"docking_tmp", "docked_ligands.mol2"))
        scorer = self.get_scorer_result()
        hs = scorer.get_hotspot()

        scored_ligs = []
        for mol in docked_mols:
            scored_l = hs.score(mol)
            ligand_score = np.mean([a.partial_charge for a in scored_l.heavy_atoms])
            scored_l.identifier += "_{}".format(round(ligand_score, 2))
            scored_ligs.append(scored_l)

        with MoleculeWriter(os.path.join(self.hotspot_path, "scored_docks.mol2")) as writer:
            for ligand in scored_ligs:
                writer.write(ligand)

    def scaled_score_ligands(self, tolerance):
        """
        Applies linear scaling to scores assigned to atom, depending on distance between atom and the scored point.
        :param int tolerance: How many gridpoints away is it acceptable for an atom to be from the nearest point of its corresponding map.
        :return: 
        """
        dsc = self.get_scorer_result()
        hs = dsc.get_hotspot()
        all_ligs = MoleculeReader(join(self.hotspot_path, "docking_tmp", "docked_ligands.mol2"))

        scored_ligs = []

        for lig in all_ligs:
            scored_lig = dsc.get_scaled_score(lig, tolerance, hs)[0]
            ligand_score = np.mean([a.partial_charge for a in scored_lig.heavy_atoms])
            scored_lig.identifier += "_{}".format(round(ligand_score, 2))
            scored_ligs.append(scored_lig)


        with MoleculeWriter(os.path.join(self.hotspot_path, "scored_docks.mol2")) as writer:
            for ligand in scored_ligs:
                writer.write(ligand)

if __name__ == "__main__":
    fu_path = join(os.getcwd(), "NUDT5_siteB.csv")
    hs_path = "/home/jin76872/Desktop/Mih/Data/NUDT5A/NUDT5A_fragment_hits/hotspot_results/NUDT5_0122_AB_bound_results"
    fu = FollowUp(fu_path, hs_path)
    #fu.score_docked_ligands()
    fu.scaled_score_ligands(tolerance=2)
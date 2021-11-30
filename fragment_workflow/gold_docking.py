"""
Classes that will do a docking run using GOLD via the CSD Python API
Documentation for this can be found here:
https://downloads.ccdc.cam.ac.uk/documentation/API/descriptive_docs/docking.html

Version 2.3.0 of the CSD Python API used here
"""
# Generic Python imports
from pathlib import Path
import pickle

# CCDC imports
from ccdc.docking import Docker
from ccdc.protein import Protein
from ccdc.molecule import Molecule
from ccdc.io import MoleculeWriter, MoleculeReader, _CSDDatabaseLocator
from ccdc.descriptors import MolecularDescriptors
from ccdc.conformer import MoleculeMinimiser
from ccdc.utilities import _private_importer
with _private_importer():
    import ChemicalAnalysisLib
    import ConformerGeneratorLib


class GOLD_docker():
    """
    Handles the docking of small molecules into a receptor using the CSD Python API
    """
    def __init__(self, protein_path, ligand_path, output_dir, reference_ligand_path=None,  autoscale=10.0, fitness_function='plp', minimise_ligand=False, prepare_protein=True, prepare_ligand=True):
        """
        
        :param protein_path:target protein. If a list of proteins is supplied, GOLD will perform ensemble docking
        :param ligand_path: the ligand to be docked 
        :param prepare_protein bool: Default assumes that no prior processing has been done on the structure and will protonate using the CCDC API
                                    Set to False if the protein has been previously prepared for the docking run
        :param prepare_ligand bool: Default assumes that no prior processing has been done on the ligand 
        """

        self.input_protein_path = protein_path
        self.input_ligand_path = ligand_path
        self.prepare_protein = prepare_protein
        self.prepare_ligand = prepare_ligand
        self.autoscale=autoscale
        self.minimise_ligand = minimise_ligand
        self.fitness_function = fitness_function
        self.results_directory =  output_dir

        if not Path(self.results_directory).exists():
            Path(self.results_directory).mkdir()

        # Assume that if no reference lligand is supplied, we are doing native docking
        if not reference_ligand_path:
            self.reference_ligand_path = ligand_path
        else:
            self.reference_ligand_path = reference_ligand_path

        self.lig_name = Path(ligand_path).name.split('.')[0]
        self.prot_name = Path(protein_path).name.split('.')[0]

        self.gold_results_directory = str(Path(output_dir, 'GOLD_docking_{}'.format(self.lig_name)).resolve())
        self.conf_file_location = str(Path(self.gold_results_directory, 'api_gold.conf').resolve())

        if not Path(self.gold_results_directory).exists():
            Path(self.gold_results_directory).mkdir()

        if self.prepare_protein:
            self.prepared_protein_path = str(Path(output_dir, 'prepared_protein_{}.mol2'.format(self.prot_name)).resolve())
        else:
            self.prepared_protein_path = self.input_protein_path

        if self.prepare_ligand:
            self.prepared_ligand_path = str(Path(output_dir, 'prepared_ligand_{}.mol2'.format(self.lig_name)).resolve())
        else:
            self.prepared_ligand_path = self.input_ligand_path


        self.docking_result=None

    def prepare_protein_for_dock(self):
        """
        
        :return: 
        """
        prot = Protein.from_file(self.input_protein_path)
        prot.identifier = self.prot_name
        prot.remove_all_waters()
        prot.remove_all_metals()
        prot.add_hydrogens()

        prot.detect_ligand_bonds()

        for l in prot.ligands:
            print(l.identifier)
            prot.remove_ligand(l.identifier)
        print('Ligands reminaing {}'.format(len(prot.ligands)))

        # Save the protein
        protwr = MoleculeWriter(self.prepared_protein_path)
        protwr.write(prot)

    @staticmethod
    def ccdc_mol_from_smiles(smiles, identifier=None, generate_initial_sites=True):
        """
        Pete's function for making a ccdc molecule with initial coordinates from a smiles string.
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

    def prepare_ligand_for_dock(self):
        """
        
        :return: 
        """
        # TODO: behaviour in case there's several ligands in the file?

        lig = MoleculeReader(self.input_ligand_path)[0]
        # Apply to our supplied ligand the same functions that ligand_prep would to a CSD entry.
        lig.identifier = self.lig_name # Note -> self.lig_name should be the name of the query ligand, not the reference (unless they are same)

        lig.remove_unknown_atoms()
        #lig.assign_bond_types()



        # Does it matter what oder you protonate and assign hydrogens in?
        Docker.LigandPreparation()._protonation_rules.apply_rules(lig._molecule)
        lig.add_hydrogens()

        # Standrdises to CSD conventions - not entirely sure if this is necessary.
        lig.standardise_aromatic_bonds()
        lig.standardise_delocalised_bonds()

        if self.minimise_ligand:
            # If the ligand has no 3D coordinates, the minimisation won't work. So let's generate some:
            if not lig.is_3d:
                print(f'Input ligand {lig.identifier} has no 3D coords. Generating 3D coords')
                lig = self.ccdc_mol_from_smiles(smiles=lig.smiles, identifier=lig.identifier)

            # Minimise the ligand conformation
            molminimiser = MoleculeMinimiser()
            lig = molminimiser.minimise(lig)

        print('Checking if ligand sucessfully minimised', type(lig))

        # Save the prepped ligand:
        ligwr = MoleculeWriter(self.prepared_ligand_path)
        ligwr.write(lig)

    def dock(self, number_poses=50):
        """
        
        :return: 
        """
        # Set up protein and ligand, in case they need to be

        if self.prepare_protein:
            self.prepare_protein_for_dock()

        if self.prepare_ligand:
            self.prepare_ligand_for_dock()

        reference_ligand = MoleculeReader(self.reference_ligand_path)[0]
        prepared_protein = Protein.from_file(self.prepared_protein_path)
        prepared_ligand = MoleculeReader(self.prepared_ligand_path)[0]

        # Set up the docking run
        docker = Docker()
        docker._conf_file_name = self.conf_file_location
        docker_settings = docker.settings
        # Prevent it from generating a ton of output ligand files - the ranked docks are in 'concat_ranked_docked_ligands.mol2'
        docker_settings._settings.set_delete_rank_files(True)
        docker_settings._settings.set_delete_empty_directories(True)
        docker_settings._settings.set_delete_all_initialised_ligands(True)
        docker_settings._settings.set_delete_all_solutions(True)
        docker_settings._settings.set_delete_redundant_log_files(True)

        # Set up the binding site. Since the sites are based on ragment hits, generate a binding site around the starting hit.
        docker_settings.reference_ligand_file = self.reference_ligand_path
        docker_settings.binding_site = docker_settings.BindingSiteFromLigand(prepared_protein,
                                                                             reference_ligand,
                                                                             6.0)
        # Default distance around ligand is 6 A. Should be ok for small fragments.Also, ALAS2 site is tiny.

        docker_settings.add_protein_file(self.prepared_protein_path)

        # Choose the fitness function: options: ['goldscore', 'chemscore', 'asp', 'plp']. plp is the default.
        docker_settings.fitness_function = 'plp'
        docker_settings.autoscale = self.autoscale
        docker_settings.early_termination = False
        docker_settings.output_directory = self.gold_results_directory
        docker_settings.output_file = str(Path(self.results_directory, 'concat_ranked_docked_ligands.mol2').resolve())

        # Add the ligand
        docker_settings.add_ligand_file(self.prepared_ligand_path,
                                        number_poses)  # Second argument determines how many poses are saved

        # Perform the docking:
        gold_result = docker.dock(file_name=self.conf_file_location)
        # pickle.dump(obj=gold_result,file=Path(self.results_directory, 'gold_result').open())
        self.docking_result = gold_result

        return gold_result

    @staticmethod
    def match_heavy_atoms(mol1, mol2):
        """
        Don't think this would work for molecules that are not identical ot even indexed identically,
        but should work for testing.
        """
        heavy1 = mol1.heavy_atoms
        heavy2 = mol2.heavy_atoms
        common = set([a.label for a in heavy1]).intersection(set([b.label for b in heavy2]))
        # print(list(common))

        pairs = []
        for c in common:
            h1 = [a for a in heavy1 if a.label == c][0]
            h2 = [b for b in heavy2 if b.label == c][0]
            pairs.append((h1, h2))
        return pairs


    def docks_to_ref_rmsd(self):
        # Only calculate for complete docking results!
        docks = [l.molecule for l in self.docking_result.ligands]
        ref_lig =  MoleculeReader(self.prepared_ligand_path)[0]
        rmsds = [MolecularDescriptors.rmsd(ref_lig, nd,
                                           exclude_hydrogens=True,
                                           atoms=self.match_heavy_atoms(ref_lig, nd)) for nd in docks]
        return rmsds

    @staticmethod
    def read_docked_ligands(docked_ligs_path):
        """
        Assumes that the conf file is in the same directory as the docked ligands. 
        Which it should be!!!
        :param docked_ligs_path: str/path.
        :return: 
        """
        conf_settings = Docker.Settings.from_file(str(Path(docked_ligs_path.parent, 'api_gold.conf').resolve()))
        docked_ligs = Docker.Results.DockedLigandReader(str(docked_ligs_path.resolve()), settings=conf_settings)

        return docked_ligs

if  __name__ == "__main__":
    ref_lig = '/dls/science/groups/i04-1/software/mihaela/Data/ALAS2A/gold_fragment_redocks/1317/A:LIG901.mol2'
    in_lig = '/dls/science/groups/i04-1/software/mihaela/Data/ALAS2A/graph_newtork_followups/x1317-suggestions/fu_0/fu_0.sdf'
    ref_rec = '/dls/science/groups/i04-1/software/mihaela/Data/ALAS2A/gold_fragment_redocks/1317/ALAS2A-x1317_1.pdb'
    test_dir = '/dls/science/groups/i04-1/software/mihaela/Scripts/repos_scripts/fragment_workflow/test_docking'

    res = GOLD_docker(protein_path=ref_rec,
                      ligand_path=in_lig,
                      output_dir=test_dir,
                      reference_ligand_path=ref_lig,
                      minimise_ligand=True)
    res.dock(number_poses=100)




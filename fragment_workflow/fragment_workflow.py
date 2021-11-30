"""
Script providing ???
"""
# Generic Python imports
from pathlib import Path
from multiprocessing import Pool

#Scientific Python imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ccdc imports
# from ccdc.io import MoleculeWriter, MoleculeReader, _CSDDatabaseLocator
# from ccdc.molecule import Molecule
# from ccdc.utilities import _private_importer
# with _private_importer():
#    import ChemicalAnalysisLib
#    import ConformerGeneratorLib

# from gold_docking import GOLD_docker


# Query the graph network
from fragalysis_preproc.data import *
from fragalysis_preproc.graph import *
import sys

def query_graph(input_smiles):
    """
    Queries the graph network for suggestions based on the smiles string provided in input_smiles
    
    :param input_smiles: the SMILES string of the query fragment
    :type input_smiles str:
    
    :return: list of smiles strings for the suggested follow-up compounds
    """
    # Initiate the object that queries the graph network
    graph_query = GraphRequest()

    # Give it the Smiles string of the query fragment
    graph_query.set_smiles_url(smiles=input_smiles)

    # Fetches the result as a json
    res_json = graph_query.get_graph_json()
    print(f"Size of Fragalysis json: {sys.getsizeof(res_json)}")

    # Flatten the json result to get the smiles strings
    flat_res = flatten_json(res_json)

    # Ask Rachael!!! Based on the key patterns, keys that end in 'end' correspond to the smiles follow-ups
    valid_smiles = []
    for k, v in flat_res.items():
        if '_end' in k:
            valid_smiles.append(v)

    return valid_smiles

def ccdc_mol_from_smiles(smiles, identifier=None, generate_initial_sites=True):
    """
    IS currently a duplicate function here and in gold docking - come up with helper module for these
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


def save_followups(valid_smiles, fu_dir, frag_name=None):
    """
    Given a list of smiles strings and a directory, saves them in individual folders.
    :param valid_smiles: 
    :param fu_dir: path to where the fragemnts get saved.
    :return: list of paths (str) to where the ligands can be found.
    """
    # Convert to a Path if not already one
    if type(fu_dir) is str:
        fu_dir = Path(fu_dir)

    # Make sure path is there
    if not fu_dir.exists(): fu_dir.mkdir()

    # List to store the paths to the saved ligands
    saved_followup_paths = []

    for i, smi in enumerate(valid_smiles):
        # Generate a ccdc molecule with 3D coordinates

        if frag_name:
            mol_name = f'{frag_name}_fu_{i}'
        else:
            mol_name = f'fu_{i}'

        fu_mol_path = Path(fu_dir, mol_name)
        if not fu_mol_path.exists(): fu_mol_path.mkdir()

        fu_mol = ccdc_mol_from_smiles(smiles=smi, identifier=mol_name, generate_initial_sites=True)
        # save the molecule:
        fu_saved = str(Path(fu_mol_path, f'{mol_name}.mol2'))

        fu_wr = MoleculeWriter(fu_saved)
        fu_wr.write(fu_mol)
        saved_followup_paths.append(fu_saved)

    return saved_followup_paths

def wrapped_gold_docking(input_dict):
   """
   
   :param input_dict: 
   :return: 
   """

   gd = GOLD_docker(protein_path=input_dict['protein_path'],
                    ligand_path=input_dict['ligand_path'],
                    reference_ligand_path=input_dict['reference_ligand_path'],
                    output_dir=input_dict['output_dir'],
                    minimise_ligand=True,
                    prepare_protein=False)

   gd.dock(number_poses=100)
   res_file = Path(gd.docking_result.settings.output_file)

   if res_file.exists():
       return {Path(input_dict['ligand_path']): res_file}

# def calculate_ligand_scores(docked_data, hotspot_path):
#     """
#     :param docked_data: list of tuples (path_to_followup: docked_ranked_followups)
#     :param hotspot_path:
#     :return:
#     """
#     docked_ligands =



def main(fragment, receptor, save_dir):
    """
    
    :param fragment: path to the fragment 
    :param receptor: path to receptor
    :param save_dir: where to save output
    :return: 
    """
    # Read in the fragment and convert to smiles:
    frag = MoleculeReader(fragment)[0]
    query_smiles = frag.smiles
    print('Read in fragment')

    # Query the graph network for follow-ups
    fu_smiles = query_graph(query_smiles)
    print('Queried graph')

    fr_name = Path(fragment).name.split('.')[0]
    docking_ligands = save_followups(fu_smiles, save_dir, fr_name)
    print('Saved ligands. Started docking')

    docking_input = [{'protein_path': receptor,
                      'ligand_path': in_path,
                      'reference_ligand_path': fragment,
                      'output_dir': str(Path(in_path).parent.resolve())} for in_path in docking_ligands]

    pool = Pool(processes=20)
    docked_results = pool.map(wrapped_gold_docking, docking_input)

    return docked_results


if __name__ == "__main__":
#     tar_frag = '/home/jin76872/Desktop/Mih/Data/ALAS2A/gold_fragment_redocks/1097/x1097_A:LIG901.mol2'
#     tar_receptor = '/home/jin76872/Desktop/Mih/Data/ALAS2A/gold_fragment_redocks/1097/prepared_protein_ALAS2A-x1097_1.mol2'
#     save_dir = '/home/jin76872/Desktop/Mih/Data/ALAS2A/graph_network_followups/x1097_ALIG901_suggestions'

    import pandas as pd

    parent_frag_smiles = "CC=1C=CC(NC=2N=CN=C3NN=CC23)=CC1"
    parent_frag_name = "ACVR1A_x1344_0B"
    parent_site = "Allosteric site"

    followup_df = pd.DataFrame(
        columns=['followup_id', 'followup_smiles', 'parent_name', 'parent_smiles', 'parent_site'])

    fu_ids = []
    fu_smiles = []
    parent_names = []
    parent_smiles = []
    sites = []

    graph_smiles = query_graph(parent_frag_smiles)
    print(len(graph_smiles))

    for it, fusm in enumerate(graph_smiles):
        fu_id = f"{parent_frag_name}_followup_{it}"

        fu_ids.append(fu_id)
        fu_smiles.append(fusm)
        parent_names.append(parent_frag_name)
        parent_smiles.append(parent_frag_smiles)
        sites.append(parent_site)

    followup_df['followup_id'] = fu_ids
    followup_df['followup_smiles'] = fu_smiles
    followup_df['parent_name'] = parent_names
    followup_df['parent_smiles'] = parent_smiles
    followup_df['parent_site'] = sites

    followup_df.to_csv('/home/jin76872/Desktop/Mih/Data/ACVR1A/graph_network_followups.csv')







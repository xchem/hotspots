"""
SCript to calculate ensembles and ligand-based grids from an ensmebles of structures from KLIFS or elsewhere
"""
# Generic Python imports
from pathlib import Path
import numpy as np
import collections
import xml.etree.ElementTree as ET

# CCDC imports
from ccdc.protein import Protein
from ccdc.molecule import Molecule
from ccdc.io import MoleculeWriter, MoleculeReader

# Hotspot imports
from hotspots.hs_io import HotspotReader, HotspotWriter
from hotspots.hs_pharmacophore import PharmacophoreModel, _Ligand
from hotspots.grid_extension import Grid, _GridEnsemble
from hotspots.hs_ensembles import EnsembleResult
from hotspots.data  import common_solvents

# RDKit imports
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold as Murcko
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity, BulkTanimotoSimilarity


def get_ligand_hs_grids(ligands):
    blank_grd = Grid.initalise_grid([a.coordinates for l in ligands for a in l.atoms])
    feature_dic = {"apolar": blank_grd.copy(),
                   "acceptor": blank_grd.copy(),
                   "donor": blank_grd.copy()}
    for lig in ligands:
        atoms = lig.heavy_atoms
        for a in atoms:
            if a.is_donor and a.is_acceptor:
                feature_dic['acceptor'].set_sphere(point=a.coordinates, radius=1, value=1, scaling='None')
                feature_dic['donor'].set_sphere(point=a.coordinates, radius=1, value=1, scaling='None')

            elif a.is_acceptor:
                feature_dic['acceptor'].set_sphere(point=a.coordinates, radius=1, value=1, scaling='None')

            elif a.is_donor:
                feature_dic['donor'].set_sphere(point=a.coordinates, radius=1, value=1, scaling='None')

            else:
                # Note that right now, all non-donors and acceptors are being labelled as apolar. Problematic?
                feature_dic['apolar'].set_sphere(point=a.coordinates, radius=1, value=1, scaling='None')
    return feature_dic


def cluster_ligands_TSNE(ligands, target_name):

    phligs = []

    for l in ligands:
        ccdc_mol = l
        rd_mol = Chem.MolFromMol2Block(l.to_string())
        if rd_mol is None:
            continue
        AllChem.Compute2DCoords(rd_mol)
        fp = AllChem.GetMorganFingerprintAsBitVect(rd_mol, 2)
        chem_id = l.identifier
        phligs.append(_Ligand(ccdc_mol, rd_mol, fp, chem_id))

    cdic = PharmacophoreModel._cluster_ligands(phligs, target_name)

    return cdic

def draw_clustered_ligands(cluster_dict, out_dir=Path('.')):
    """
    Provides an image of the ligands in each cluster.
    :param cluster_dict: output of cluster_ligands(); dictionary of {cluster:[ligands]}
    :return: 
    """
    from rdkit.Chem import Draw

    for clust, ligs in cluster_dict.items():
        # Get 2D coordinates for the ligands
        rd_ligs = [l.rdmol for l in ligs]
        for r in rd_ligs: AllChem.Compute2DCoords(r)

        img = Draw.MolsToGridImage(rd_ligs, legends= [l.chemical_id for l in ligs])
        img.save(Path(out_dir, f'cluster_{clust}.png'))

def ligand_grids_from_all_ligands(ligands, target_name, cluster_mode):
    """
    
    :param ligands: ccdc molecules
    :param target_name: str, name of the ensmeble the ligands are for
    :param cluster_mode: whether to use Pete's clustering method ('TSNE'), or Tanimoto ('Tanimoto') 
    :return: dictionary of ligand grids.
    """

    if cluster_mode == 'TSNE':
        cluster_dict = cluster_ligands_TSNE(ligands, "{}_{}".format(target_name, cluster_mode))

    else:
        # TODO: Put the Tanimoto clustering here!
        cluster_dict = {}

    # Generate the ligand amps and save them
    grid_path = Path(".", f"{target_name}_{cluster_mode}_ligand_grids")
    if not grid_path.exists(): grid_path.mkdir()

    image_dir = Path(grid_path, "clustered_ligands")
    if not image_dir.exists(): image_dir.mkdir()

    draw_clustered_ligands(cluster_dict=cluster_dict, out_dir=image_dir)

    # Get the representative ligands for the clusters - for now, just take the first ligand
    reps_ccdc = [l[0].ccdc_mol for l in cluster_dict.values() if len(l) != 0]
    rep_lig_path = Path(grid_path, "representative_ligands.sdf")

    with MoleculeWriter(str(rep_lig_path.resolve())) as wr:
        for r in reps_ccdc:
            wr.write(r)

    lig_grids = get_ligand_hs_grids(reps_ccdc)
    for key, val in lig_grids.items():
        val.write(str(Path(grid_path, f"ligands_{target_name}_{key}.ccp4").resolve()))

    return reps_ccdc



def get_grid_from_plip_coords(probe, coords, out_path, padding=4.0):
    """
    
    :param probe: 
    :param coords: 
    :param padding: 
    :return: 
    """
    Coordinates = collections.namedtuple('Coordinates', ['x', 'y', 'z'])

    dims = np.array(coords)
    print(dims.shape)
    min_coords = np.array([np.min(dims[:, 0]), np.min(dims[:, 1]), np.min(dims[:, 2])])
    max_coords = np.array([np.max(dims[:, 0]), np.max(dims[:, 1]), np.max(dims[:, 2])])

    origin = Coordinates(x=round(min_coords[0] - padding),
                         y=round(min_coords[1] - padding),
                         z=round(min_coords[2] - padding))

    far_corner = Coordinates(x=round(max_coords[0] + padding),
                             y=round(max_coords[1] + padding),
                             z=round(max_coords[2] + padding))

    plip_grid = Grid(origin=origin, far_corner=far_corner, spacing=0.5, default=0, _grid=None)

    for coo in coords:
        plip_grid.set_sphere(point=coo, radius=1, value=1, scaling='None')

    plip_grid.write(str(Path(out_path, f"plip_{probe}.ccp4").resolve()))

    return plip_grid


def get_plip_coordinates(report_path, ligand_id=None):
    """
    
    :param report_path: Path to the output of PLIP (in xml format)
    :param ligand_id: 
    :return: dictinary of {'probe': coordinates}
    """
    tree = ET.parse(report_path)
    root = tree.getroot()

    hydrophobes = []
    lig_acceptors = []
    lig_donors = []

    solvents = common_solvents()

    # There might be multiple binding sites. Need to make sure we're only getting coordinates out for the main binding site.
    for binding_site in root.findall('bindingsite'):

        # Go away glycerol!
        nest = binding_site.find('identifiers')
        lig_code = nest.find('hetid')

        #if str(lig_code.text) == str(ligand_id):
        if str(lig_code.text).upper() not in solvents:

            # Find all the interactions
            interactions = binding_site.find('interactions')

            # Extract the hydrophobic ones
            hydrophobic_inters = interactions.find('hydrophobic_interactions')

            for h in hydrophobic_inters:
                lig_coords = [float(x.text) for x in h.find('ligcoo')]
                hydrophobes.append(lig_coords)

            # Get interactions labelled as H-bonds
            h_bonds = interactions.find('hydrogen_bonds')

            for hb in h_bonds:
                print(hb.tag, hb.attrib)
                # Find the ligand coordinates
                lig_coords = [float(x.text) for x in hb.find('ligcoo')]

                # check if the protein is acting as accetor or donor. Ligand should be the opposite.
                protisdon = hb.find('protisdon').text
                print(protisdon, type(protisdon))

                if protisdon == "True":
                    # ... then ligisacc
                    lig_acceptors.append(lig_coords)
                else:
                    lig_donors.append(lig_coords)

    return {'apolar': hydrophobes,
            'acceptor': lig_acceptors,
            'donor': lig_donors}

def make_plip_grids(report_paths, out_path):
    """
    
    :param report_paths_lig_dict: dictionary of form {Path(report.xml}:lig_id}
    :param out_path: where to save the PLIP grids
    :return: grid_list of form {'probe':grid} -> can be fed into a hotspot result
    """

    # Store all the PLIP coordinates accross all of the ligand-bound structures for a protein
    master_coords_dict = {}

    for rep_path in report_paths:
        struct_PLIP_dict = get_plip_coordinates(report_path=rep_path)

        for probe, coords_list in struct_PLIP_dict.items():
            try:
                master_coords_dict[probe].extend(coords_list)

            except KeyError:
                master_coords_dict[probe] = coords_list

    # Make grids from the coordinates
    plip_grid_dict = {}

    for p, coords in master_coords_dict.items():
        plip_grid_dict[p] = get_grid_from_plip_coords(p, coords, out_path)

    return plip_grid_dict



if __name__ == "__main__":

    import sys

    sys.path.insert(0, "/home/jin76872/Desktop/Mih/Scripts/repos_scripts/ensemble_pipeline/")
    from run_hotspots import ParalleliselRunner
    import luigi
    import tempfile
    import shutil

    tempfile.tempdir = '/home/jin76872/Desktop/Mih/Data/tmp_superstar_ghecom'



    ensemble_path = Path("/home/jin76872/Desktop/Mih/Data/ensemble_maps_validation/p38a/p38a")
    # ens_ligands = [MoleculeReader(str(p.resolve()))[0] for p in ensemble_path.glob("*/ligand.mol2")]
    #
    # repr_ligs = ligand_grids_from_all_ligands(ligands=ens_ligands,
    #                                           target_name='p38a',
    #                                           cluster_mode='TSNE')
    #
    # # Get the PDB IDs of the representative liagnds:
    # repr_pdbs = [r.identifier for r in repr_ligs]
    #
    # pdb_paths = [p for p in ensemble_path.glob("*/protein.mol2") for r in repr_pdbs if f'{r.split(".")[0]}' in str(p.resolve()) and f'_chain{r.split(".")[1].upper()}'in str(p.resolve())]
    #
    # # Get rid of all the alts unless they are AltA. Not sure if there is an altC anywhere, but just in case
    #
    # pdb_paths = [p for p in pdb_paths if ('altB' or 'altC') not in str(p.resolve())]
    #
    # # Copy the protein paths and save them
    # reduced_dir = Path(ensemble_path, "reduced_ensemble")
    # if not reduced_dir.exists(): reduced_dir.mkdir()
    #
    # prot_paths = []
    # for p in pdb_paths:
    #     cop_path = Path(reduced_dir, p.parent.name)
    #     if not cop_path.exists(): cop_path.mkdir()
    #     prot_path = Path(cop_path, f'{p.parent.name}.mol2')
    #     if not prot_path.exists():
    #         #prot_path.symlink_to(p)
    #         shutil.copy(p, prot_path)
    #     prot_paths.append(str(prot_path.resolve()))
    #
    # luigi.build([ParalleliselRunner(in_list = prot_paths,
    #                                nrot=3000,
    #                                charged=False,
    #                                data_source='KLIFS')], local_scheduler=True, workers=20)

    hs_results_paths = ensemble_path.glob("*/fullsize_hotspots_3000/out.zip")
    repr_ligands = MoleculeReader('/home/jin76872/Desktop/Mih/Data/ensemble_maps_validation/p38a_TSNE_ligand_grids/representative_ligands.sdf')

    shrunk_results = shrink_hotspot_maps(hs_result_paths=hs_results_paths, ligands=repr_ligands, padding=4.0)

    # Create the ensmeble maps:
    # Set up parameters for the Ensemble Maps. Shown below are the default values, but the ensemble_settings allow for tweaking them.

    ensemble_settings = EnsembleResult.Settings()
    # Take the median value across the ensemble at each point in space
    ensemble_settings.combine_mode = 'median'

    # Include all points in the apolar maps
    ensemble_settings.apolar_frequency_threshold = 0.0

    # For the polar maps, use only points with nonzero scores in over 20% of the ensemble.
    ensemble_settings.polar_frequency_threshold = 20.0









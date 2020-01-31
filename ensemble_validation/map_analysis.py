#Generic imports
from pathlib import Path
import collections
import xml.etree.ElementTree as ET

# Scientific Python imports
import numpy as np
from sklearn.metrics import precision_recall_curve, average_precision_score, auc
import pandas as pd
import matplotlib.pyplot as plt

# Hotspot imports
from hotspots.grid_extension import Grid, _GridEnsemble
from hotspots.hs_io import HotspotReader, HotspotWriter
from hotspots.hs_ensembles import EnsembleResult
from hotspots.data  import common_solvents


def get_grid_from_plip_coords(probe, coords, out_path, padding=4.0):
    """

    :param probe: 
    :param coords: 
    :param padding: 
    :return: 
    """
    Coordinates = collections.namedtuple('Coordinates', ['x', 'y', 'z'])

    dims = np.array(coords)
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

        # if str(lig_code.text) == str(ligand_id):
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
                # Find the ligand coordinates
                lig_coords = [float(x.text) for x in hb.find('ligcoo')]

                # check if the protein is acting as accetor or donor. Ligand should be the opposite.
                protisdon = hb.find('protisdon').text

                if protisdon == "True":
                    # ... then ligisacc
                    lig_acceptors.append(lig_coords)
                else:
                    lig_donors.append(lig_coords)

    return {'apolar': hydrophobes,
            'acceptor': lig_acceptors,
            'donor': lig_donors}


def make_plip_grids(report_paths, out_path, padding=6.0):
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
        plip_grid_dict[p] = get_grid_from_plip_coords(p, coords, out_path, padding)

    return plip_grid_dict

def calculate_precision_recall(query_grid, reference_grid, threshold=0):
    """
    
    :param query_grid: Hotspots.grid_extension.Grid
    :param reference_grid: Hotspots.grid_extension.Grid
    :return: 
    """
    reference_mask = reference_grid > threshold

    query_vals = np.array(query_grid.to_vector())
    true_labels = np.array(reference_mask.to_vector())
    precision, recall, thresholds = precision_recall_curve(true_labels, query_vals)

    return{'precision':precision,
           'recall': recall,
           'thresholds':thresholds}

def plot_precision_recall(query_grid, reference_grid, probe_name):
    """
    
    :param query_grid: 
    :param reference_grid: 
    :return: 
    """
    res_dict = calculate_precision_recall(query_grid, reference_grid)
    plt.plot(res_dict['recall'], res_dict['precision'])
    auc_ = auc(res_dict['recall'], res_dict['precision'])

    # Add AUC and probe information to the plot
    props = dict(boxstyle='round', alpha=0.5)
    ax = plt.gca()
    ax.text(0.05, 0.95, f'Probe: {probe_name}, \n AUC: {round(auc_, 2)}', transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)


def plot_precision_recall_figure(query_grids, reference_grids):
    """
    
    :param query_grids: 
    :param reference_grids: 
    :return: 
    """
    probes = ['donor', 'acceptor', 'apolar']

    fig1=plt.figure(figsize=(6, 8))
    for n, key in enumerate(probes):
        plt.subplot(3, 1, n + 1)
        plot_precision_recall(query_grids[key], reference_grids[key], key)
        if n < 2:
            plt.xticks([])
        if n == 2:
            plt.xlabel("Recall")
        if n == 1:
            plt.ylabel("Precision")
    #plt.savefig(join(stem, "{}_mean_spread.png".format(prot_name)))
    plt.show()
    #plt.close()


def shrink_hotspot_maps(hs_result_paths, ligands, padding=4.0):
    """
    Given the list of ligands in the ensemble and some hotspot maps, will shrink all the maps in the 
    :param hs_result_paths: a list of Paths to precalculated hotspot results. Should be all for the same target (or actually, targets that we are looking to compare. 
    :param ligands: a list of ccdc molecules corersponding to the ensmeble ligands. Needed to define the binding site of interest.
    :return: a list of *shrunk* hotspot results
    """

    # Find the largest ligand and use it to define the binding site
    mws = [l.molecular_weight for l in ligands]
    biggest_lig = ligands[mws.index(max(mws))]

    # Get the dimensions in space of the largest ligand
    dims = np.array([a.coordinates for a in biggest_lig.atoms])
    min_coords = np.array([np.min(dims[:, 0]), np.min(dims[:, 1]), np.min(dims[:, 2])])
    max_coords = np.array([np.max(dims[:, 0]), np.max(dims[:, 1]), np.max(dims[:, 2])])

    # Add some padding in both directions:
    min_coords -= padding
    max_coords += padding

    # Now shrink all the hotspot grids to the min and max dimensions
    shrunk_hs_results = []

    for hpath in hs_result_paths:
        hs_res = HotspotReader(str(hpath.resolve())).read()
        probes = hs_res.super_grids.keys()

        # now to shrink the grids for each probe
        for p in probes:
            hs_res.super_grids[p] = EnsembleResult.shrink_to_binding_site(in_grid=hs_res.super_grids[p],
                                                                          new_origin=min_coords,
                                                                          new_far_corner=max_coords)
            shrunk_hs_results.append(hs_res)

            h_out_dir = Path(hpath.parent, 'binding_site_maps')
            if not h_out_dir.exists(): h_out_dir.mkdir()
            with HotspotWriter(str(h_out_dir.resolve()), visualisation="pymol", grid_extension=".ccp4", zip_results=False) as writer:
                writer.write(hs_res)

    return shrunk_hs_results




"""
Miscellaneous utility functions for making and handling hotspots maps.

Main structures:
struct_dict = {['prot_name']:(path_to_protein, ccdc.protein.Protein object)}
grid_dict = {['probe']: ccdc.utilities.Grid object}
"""

from __future__ import print_function, division
from fragment_hotspot_maps.fragment_hotspot_maps import Hotspots, HotspotsHelper
from ccdc.protein import Protein
from ccdc.utilities import Grid
import os
from os.path import join
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import atexit


# Functions for reading in grids and proteins

def read_proteins(stem, substring):
    """
    Reads all pdb files in a directory containing the substring
    :param stem: str(path to directory)
    :param substring: str(name of protein)
    :return: Python dictionary object
    """
    struct_dict = {}
    for root, subdirs, files in os.walk(stem):
        for filename in files:
            if substring in filename:
                if '.pdb' in filename:
                    name_path = os.path.join(root, filename)
                    prot_name = filename.replace('.pdb', '')
                    struct_dict[prot_name] = (name_path, Protein.from_file(name_path))
    return struct_dict


def read_grids(stem, substring, charged_probes=False):
    """
    Creates a grid_dic from a given directory
    :param stem: str(path to directory)
    :param substring: str(protein name)
    :param charged_probes: whether o look for positive and negative grids
    :return: Python dictionary
    """
    grid_dic = {}
    if charged_probes == False:
        probe_list = ['apolar', 'acceptor', 'donor']
    else:
        probe_list = ['apolar', 'acceptor', 'donor', 'positive', 'negative']

    for root, subdirs, files in os.walk(stem):
        for filename in files:
            if substring in filename:
                name_path = os.path.join(root, filename)
                if '.ccp4' in filename:
                    prot_name = filename.replace('.ccp4', '')
                elif '.grd' in filename:
                    prot_name = filename.replace('.grd', '')
                else:
                    continue
                for probe in probe_list:
                    if probe in filename:
                        grid_dic[probe] = Grid.from_file(name_path)

    return grid_dic


def find_grid(stem, prot_name, probe=None):
    """
    Loops through a directory to find a specific grid. Doesn't discriminate probes
    :param prot_name: 
    :return: ccdc.utilities.Grid object
    """
    for root, subdirs, files in os.walk(stem):
        for filename in files:
            if prot_name in filename and '.ccp4' in filename:
                name_path = os.path.join(root, filename)
                if probe != None:
                    if probe in filename:
                        grid = Grid.from_file(name_path)
                        return grid
                else:
                    grid = Grid.from_file(name_path)
                    return grid


# Functions for generating Hotspot Maps

def prepare_protein(prot_name, struct_dict, ch='A'):
    """
    Function to prepare proteins for a SuperStar Run.
    :param prot_name: struct_dict key
    :param struct_dict: 
    :param ch: str(chain). If not specified, assumes chain 'A'
    :return: Python dictionary
    """
    prot = struct_dict[prot_name][1]
    for chain in prot.chains:
        print(chain.identifier)
        if chain.identifier != ch:
            prot.remove_chain(chain.identifier)
    prot.remove_all_waters()
    for ligand in prot.ligands:
        prot.remove_ligand(ligand.identifier)
    prot.remove_all_metals()
    prot.add_hydrogens()
    return struct_dict


def make_hotspot_maps(stem, substring, grid_out_dir, sele='A', ghecom_exec=None, ite=1, compressed=False):
    """
    Creates Hotspot Maps for proteins in the stem directory
    :param stem: str(path to dir)
    :param substring: str(protein name)
    :param grid_out_dir: str(path to where maps are saved)
    :param sele: str(which chain to get maps for)
    :param ghecom_exec: str(path to Ghecom executable. Optional)
    :param ite: int(how many maps per structure)
    :param compressed : bool, whether to run compress_grid() before saving.
    :return: fragment_hotspot_maps.HotspotResults object
    """
    struct_dict = read_proteins(stem, substring)
    for prot_name in struct_dict:
        print(prot_name)
        for i in range(ite):
            prepare_protein(prot_name, struct_dict, sele)
            h = Hotspots()
            result = h.from_protein(struct_dict[prot_name][1], fname=struct_dict[prot_name][0],
                                    ghecom_executable=ghecom_exec, charged_probes=False)
            for probe, g in result.super_grids.items():
                if compressed == True:
                    g = compress_grid(g)
                g.write(join(grid_out_dir, (prot_name + ('_{}_{}.ccp4'.format(probe, str(ite))))))
            atexit._run_exitfuncs()
    return result


# Plotting Histograms

def get_grid_dic_histograms(grid_dic, out_dir, prot_name, suffix=None):
    """
    Same as histograms function in main hotspots code, but for any grid_dic
    :param grid_dic: Python dictionary {[probe]:ccdc.Utilities.Grid}
    :param out_dir: str(path to out_dir)
    :param prot_name: str
    :param suffix: str
    :return: 
    """
    data = {}
    for g in grid_dic.keys():
        grd = grid_dic[g]
        nx, ny, nz = grd.nsteps
        data[g] = np.array([grd.value(i, j, k) for i in xrange(0, nx) for j in xrange(0, ny) for k in xrange(0, nz) if
                            round(grd.value(i, j, k)) != 0])

    plt.figure(1)
    for n, key in enumerate(data.keys()):
        if suffix != None:
            np.savetxt(join(out_dir, 'data_{}_{}.txt'.format(key, suffix)), data[key])
        print(key)
        print(n)
        plt.subplot(3, 1, (n + 1))
        # hotspot_result._histogram_info(data, key, n)
        colour_dict = {"acceptor": "r", "donor": "b", "apolar": "y"}
        hist, bin_edges = np.histogram(data[key], bins=range(0, 40), normed=True)
        plt.bar(bin_edges[:-1], hist, width=1, color=colour_dict[key])
        plt.xlim(min(bin_edges), max(bin_edges))
        plt.ylim(0, 0.35)
        plt.yticks([])
        if n == 0:
            plt.title("Fragment hotspot Maps")
        if n < 2:
            plt.xticks([])
        if n == 2:
            plt.xlabel("Fragment hotspot score")
        if n == 1:
            plt.ylabel("Frequency")
    if suffix != None:
        plt.savefig(join(out_dir, (prot_name + suffix)))
    else:
        plt.savefig(join(out_dir, prot_name))
    plt.close()


def make_histogram(grd, out_dir, suffix):
    """
    Makes a histogram of a single grid
    :param grd: ccdc.utilities Grid
    :param out_dir: path (str)
    :param suffix: str
    :return: 
    """
    nx, ny, nz = grd.nsteps
    hist_arr = np.array([grd.value(i, j, k) for i in xrange(0, nx) for j in xrange(0, ny) for k in xrange(0, nz) if
                         round(grd.value(i, j, k)) > 5])
    plt.figure(1)
    plt.hist(hist_arr, bins=10)
    # plt.figtext(0.6, 0.8, ('Number of zero values:' + str(num_zeros)))
    plt.xlabel('Hotspot_score')
    plt.ylabel('Frequency')
    plt.savefig(join(out_dir, 'Hotspot_histogram_{}'.format(suffix)))
    plt.close()


# Utility functions for CCDC Grids

def grid_to_numpy(g):
    """
    Makes a ccdc.utilities Grid object into a numpy array
    :param grid:  ccdc.utilities Grid object
    :return: numpy array
    """

    def get_value(value, x, y, z):
        return value + g.value(x, y, z)

    vget_value = np.vectorize(get_value)
    nx, ny, nz = g.nsteps
    arr = np.zeros((nx, ny, nz), dtype=float)
    index_array = np.indices(arr.shape)

    arr = vget_value(arr, index_array[0], index_array[1], index_array[2])
    return arr


def are_same(g1, g2):
    """
    takes 2 grids and checks if they're different
    :param g1: ccdc.utilities Grid object
    :param g2: ccdc.utilities Grid object
    :return: 
    """

    if g1.bounding_box[0] != g2.bounding_box[0]:
        print('origins differ')

    diff_dict = [g1.nsteps[i] - g2.nsteps[i] for i in range(3)]
    for diff in diff_dict:
        if diff != 0:
            print('Dimensions are not the same')

    diff_gr = g1 - g2
    if diff_gr.count_grid() != 0:
        print('Grids not the same')
    else:
        print("They're the same")


def vrun_gaussian(g, sigma):
    """
    gaussian smoothing function, method of reducing noise in output
    Numpy vectorised version of the _run_gaussian fx in the main Hotspots Code
    :param g : (ccdc grid object)
    :param sigma: (sigma i each direction + order of the Gaussian)
    :return:new ccdc grid
    """

    def get_value(value, x, y, z):
        return value + g.value(x, y, z)

    vget_value = np.vectorize(get_value)

    nx, ny, nz = g.nsteps
    scores = np.zeros((nx, ny, nz, 1))
    index_array = np.indices(scores.shape)

    scores = vget_value(scores, index_array[0], index_array[1], index_array[2])
    mod = ndimage.filters.gaussian_filter(scores, sigma=sigma)

    def set_value(x, y, z):
        new_grid.set_value(x, y, z, mod[x, y, z, 0])

    vset_value = np.vectorize(set_value)

    # equivalent of 'copy_and_clear' from main code
    new_grid = g.copy()
    new_grid *= 0

    vset_value(index_array[0], index_array[1], index_array[2])

    return new_grid


# Saving Grids

def compress_grid(grid):
    """
    Compresses a grid.
    :param grid: ccdc.utilities Grid
    :return: ccdc.utilities Grid
    """

    arr = grid_to_numpy(grid)
    ind_arr = np.where(arr>0)
    or_diff = [min(ind_arr[i]) for i in range(3)]
    print(or_diff)
    far_diff = [grid.nsteps[i] - max(ind_arr[i]) for i in range(3)]
    print(far_diff)
    if int(min(or_diff)) < 3 or int(min(far_diff)) < 3:
        p = 0
    else:
        p = 2
    region = (min(ind_arr[0])-p, min(ind_arr[1])-p, min(ind_arr[2])-p, max(ind_arr[0])+p, max(ind_arr[1])+p, max(ind_arr[2])+p)
    compr_grid = Grid.sub_grid(grid, region)

    return compr_grid


# Functions to get selectivity maps between 2 grids

def divide_select(g1, g2):
    """
    Makes a selectivity map by dividing two grids
    :param g1: ccdc.utilities.Grid
    :param g2: ccdc.utilities.Grid
    :return: ccdc.utilities Grid
    """

    h = HotspotsHelper()
    g1, g2 = h._common_grid(g1, g2)
    com_bound_box = g1.bounding_box
    com_spacing = g1.spacing

    arr1 = grid_to_numpy(g1)
    arr1[arr1 < 1] = 1.0
    arr2 = grid_to_numpy(g2)
    arr2[arr2 < 1] = 1.0

    sel_arr = np.divide(arr1, arr2)
    sel_arr[sel_arr == 1.0] = 0.0

    sel_map = Grid(origin=com_bound_box[0], far_corner=com_bound_box[1], spacing=com_spacing)

    for x in range(sel_map.nsteps[0]):
        for y in range(sel_map.nsteps[1]):
            for z in range(sel_map.nsteps[2]):
                sel_map.set_value(x, y, z, sel_arr[x][y][z])

    return sel_map


def subtract_select(g1, g2):
    """
    Makes a selectivity map by subtracting 2 grids
    :param g1: 
    :param g2: 
    :return: 
    """
    h = HotspotsHelper()
    g1, g2 = h._common_grid(g1, g2)
    sel_map = g1 - g2
    return sel_map

def filter_map(g1, g2):
    """
    Sets 2 grids to the same size and coordinate frames. Points that are zero in one grid but sampled in the other are
    set to the mean of their nonzero neighbours.
    :param g1: ccdc.utilities.Grid
    :param g2: ccdc.utilities.Grid
    :return: ccdc.utilities.Grid
    """

    def filter_point(x, y, z):
        loc_arr = np.array([g[x+i][y+j][z+k] for i in range(-1, 2) for j in range(-1, 2) for k in range(-1, 2)])
        if loc_arr[loc_arr>0].size != 0:
            print(np.mean(loc_arr[loc_arr > 0]))
            new_grid[x][y][z] = np.mean(loc_arr[loc_arr>0])

    vfilter_point = np.vectorize(filter_point)
    h = HotspotsHelper()
    g1 = compress_grid(g1)
    g2 = compress_grid(g2)
    g1, g2 = h._common_grid(g1, g2, padding=1)
    com_bound_box = g1.bounding_box
    com_spacing = g1.spacing

    arr1 = grid_to_numpy(g1)
    arr2 = grid_to_numpy(g2)

    b_arr1 = np.copy(arr1)
    b_arr2 = np.copy(arr2)

    b_arr1[b_arr1>0] = 1.0
    b_arr2[b_arr2>0] = -1.0

    diff_arr = b_arr1 + b_arr2

    unmatch1 = np.where(diff_arr == 1)
    unmatch2 = np.where(diff_arr == -1)

    g = arr1
    new_grid = np.copy(arr1)
    vfilter_point(unmatch2[0], unmatch2[1], unmatch2[2])
    f_arr1 = np.copy(new_grid)
    f_arr1[f_arr1<1] = 1

    g = arr2
    new_grid = np.copy(arr2)
    vfilter_point(unmatch1[0], unmatch1[1], unmatch1[2])
    f_arr2 = np.copy(new_grid)
    f_arr2[f_arr2<1] = 1

    sel_arr = np.divide(f_arr1, f_arr2)
    sel_arr[sel_arr == 1] = 0

    sel_map = Grid(origin=com_bound_box[0], far_corner=com_bound_box[1], spacing=com_spacing)

    for x in range(sel_map.nsteps[0]):
        for y in range(sel_map.nsteps[1]):
            for z in range(sel_map.nsteps[2]):
                sel_map.set_value(x, y, z, sel_arr[x][y][z])
    return sel_map


from __future__ import print_function, division
from EnsemblePipeline import EnsembleIO
from GridEnsemble import GridEnsemble
from hotspots.grid_extension import Grid
from glob import glob
import os
from os.path import join, exists
import pickle
import matplotlib.pyplot as plt
import numpy as np


def create_gridensembles(ens, io):
    """
    
    :param str root_dir: Path to the pipeline root 
    :return: 
    """
    # Create and save the GridEnsembles for all of our bromodomains
    probes = ["donor", "acceptor", "apolar"]

    probe_dict = {}
    for probe in probes:
        paths = [join(io.binding_site_hotspot_paths[ens][i], "{}.ccp4".format(probe))
                 for i in range(len(io.binding_site_hotspot_paths[ens]))]
        ge = GridEnsemble(paths)
        ge.get_ensemble_array()
        save_path = join(io.ensemble_dirs[ens], "{}_{}.p".format(ens, probe))
        ge.save_gridensemble(save_path)
        probe_dict[probe] = save_path
    io.ensemble_maps[ens] = probe_dict
    io.update()

def get_p_value_maps(io):
    """
    
    :return: 
    """
    probes = ["donor", "acceptor", "apolar"]
    for probe in probes:
        ge1 = pickle.load(open(io.ensemble_maps[e1][probe], "r"))
        ge2 = pickle.load(open(io.ensemble_maps[e2][probe], "r"))

        iarr = ge1.get_difference_map(ge2, tolerance=0)
        gr = ge1.save_grid(iarr)
        gr.write(join(io.params.pipeline_root, "{}_{}_{}_ranges.ccp4".format(e1, e2, probe)))


def get_max_grids(io):
    max_dict = {}
    probes = ["donor", "acceptor", "apolar"]
    for e in io.ensembles:
        for probe in probes:
            ge = pickle.load(open(io.ensemble_maps[e][probe], "r"))
            max_g = ge.make_max_grid()
            max_dict[e] = max_g
            max_g.write(join(io.ensemble_dirs[e], "{}_{}_max.ccp4".format(e, probe)))

def get_med_grids(io):
    probes = ["donor", "acceptor", "apolar"]
    for e in io.ensembles:
        for probe in probes:
            ge = pickle.load(open(io.ensemble_maps[e][probe], "r"))
            g = ge.make_median_grid()
            g.write(join(io.ensemble_dirs[e], "{}_{}_median.ccp4".format(e, probe)))

def get_mean_grids(io):
    probes = ["donor", "acceptor", "apolar"]
    for e in io.ensembles:
        for probe in probes:
            ge = pickle.load(open(io.ensemble_maps[e][probe], "r"))
            g = ge.make_mean_grid()
            g.write(join(io.ensemble_dirs[e], "{}_{}_mean.ccp4".format(e, probe)))

def make_max_difference_maps(io):
    probes = ["donor", "acceptor", "apolar"]
    for probe in probes:
        g1 = Grid.from_file(join(io.ensemble_dirs[e1], "{}_{}_max.ccp4".format(e1, probe)))
        g2 = Grid.from_file(join(io.ensemble_dirs[e2], "{}_{}_max.ccp4".format(e2, probe)))
        diff_g = g1 - g2
        diff_g.write(join(io.params.pipeline_root, "diff_{}_{}_{}.ccp4").format(e1, e2, probe))


def make_median_difference_maps(io):
    probes = ["donor", "acceptor", "apolar"]
    for probe in probes:
        g1 = Grid.from_file(join(io.ensemble_dirs[e1], "{}_{}_median.ccp4".format(e1, probe)))
        g2 = Grid.from_file(join(io.ensemble_dirs[e2], "{}_{}_median.ccp4".format(e2, probe)))
        diff_g = g1 - g2
        diff_g.write(join(io.params.pipeline_root, "median_diff_{}_{}_{}.ccp4").format(e1, e2, probe))


def make_thresholded_maps(io, e1, e2):
    """
    
    :param io: EnsembleIO instance. 
    :param str ens1: name of reference ensemble 
    :param str ens2: name of off-target ensemble 
    :return: 
    """
    t_dir = join(io.params.pipeline_root, "thresholded_hotspot_maps")
    if not exists(t_dir):
        os.mkdir(t_dir)
    print(t_dir)
    probes = ["donor", "acceptor", "apolar"]

    for probe in probes:
        ge1 = pickle.load(open(io.ensemble_maps[e1][probe], "r"))
        ge2 = pickle.load(open(io.ensemble_maps[e2][probe], "r"))

        iarr = ge1.get_difference_map(ge2, tolerance=0)
        diff = Grid.from_file(join(io.params.pipeline_root, "diff_{}_{}_{}.ccp4").format(e1, e2, probe)).get_array()

        over3 = (iarr > 3) * diff
        gover3 = ge1.save_grid(over3)
        gover3.write(join(t_dir, "diff_{}_{}_{}_over3.ccp4".format(e1, e2, probe)))

        under3 = (iarr < 3) * diff
        gunder3 = ge1.save_grid(under3)
        gunder3.write(join(t_dir, "diff_{}_{}_{}_under3.ccp4".format(e1, e2, probe)))


def make_median_thresholded_maps(io, e1, e2, threshold):
    """

    :param io: EnsembleIO instance. 
    :param str ens1: name of reference ensemble 
    :param str ens2: name of off-target ensemble 
    :return: 
    """
    t_dir = join(io.params.pipeline_root, "median_thresholded_hotspot_maps")
    if not exists(t_dir):
        os.mkdir(t_dir)
    print(t_dir)
    probes = ["donor", "acceptor", "apolar"]

    for probe in probes:
        ge1 = pickle.load(open(io.ensemble_maps[e1][probe], "r"))
        ge2 = pickle.load(open(io.ensemble_maps[e2][probe], "r"))

        iarr = ge1.get_difference_frequency_map(other=ge2, threshold=threshold)
        gover3 = ge1.save_grid(iarr)
        gover3.write(join(t_dir, "diff_median_freq_thresh{}_{}_{}_{}.ccp4".format(threshold, e1, e2, probe)))

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
        print(key)
        print(n)
        plt.subplot(3, 1, (n + 1))
        colour_dict = {"acceptor": "r", "donor": "b", "apolar": "#FFF176"}
        plt.hist(data[key], bins=40, edgecolor='k', alpha=0.8, color=colour_dict[key])
        plt.xlim(0, 35.0)
        if n == 0:
            plt.title("Fragment hotspot Maps")
        if n == 2:
            plt.xlabel("Fragment hotspot score")
        if n == 1:
            plt.ylabel("Frequency")
    if suffix != None:
        plt.savefig(join(out_dir, "scores_{}_{}".format(prot_name, suffix)))
    else:
        plt.savefig(join(out_dir, prot_name))
    plt.close()


def plot_distribution_histograms(io, mode):
    probes = ['acceptor', 'apolar', 'donor']
    for e in io.ensembles:
        grid_dict = {}
        for probe in probes:
            paths = glob(join(io.ensemble_dirs[e], '*.ccp4'))
            print(paths)
            path = [p for p in paths if probe in p and mode in p][0]
            print(path)
            grid_dict[probe] = Grid.from_file(path)

        get_grid_dic_histograms(grid_dict, out_dir = io.ensemble_dirs[e], prot_name=e, suffix=mode)


if __name__ == "__main__":
    import json
    tar_path = "/home/jin76872/Desktop/Mih/Data/pipeline_dev/KLIFS_more"
    #create_gridensembles(tar_path)
    # KLIFS_dict_path = join(tar_path, "shrunk_hot_paths_100000.json")
    # with open(KLIFS_dict_path, "r") as f:
    #     KLIFS_dict = json.load(f)

    io = EnsembleIO(tar_path)
    ensemble_names = ['p38a', 'Erk2']
    io.ensembles = ensemble_names
    io.ensemble_dirs = {e:join(tar_path, e) for e in ensemble_names}
    #io.binding_site_hotspot_paths = KLIFS_dict

    e2 = "p38a"
    e1 = "Erk2"
    make_median_thresholded_maps(io, e1, e2, threshold=20)
    # for n in ensemble_names:
    #     create_gridensembles(n, io)
    get_p_value_maps(io)
    get_max_grids(io)
    get_med_grids(io)
    get_mean_grids(io)
    make_median_difference_maps(io)
    make_max_difference_maps(io)
    make_thresholded_maps(io, e1, e2)

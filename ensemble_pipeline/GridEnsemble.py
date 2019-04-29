from __future__ import print_function, division
from hotspots.grid_extension import Grid
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.feature_extraction.image import img_to_graph
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
import pickle
from decimal import Decimal
from scipy.stats import ks_2samp


class GridEnsemble:
    """
    Slightly different way of compiling the Ensemble
    """

    def __init__(self, path_list, grid_spacing=0.5):
        """

        :param path_list: list of paths to hotspot maps for the same probe type
        :param grid_spacing: default if 0.5, also used by Superstar. Only change if consistently working with different-spaced grids.
        """
        # Use the paths rather than loading in all grids at once - reduces memory usage for large ensembles.
        self.paths = path_list
        # Set the spacing before compiling the ensemble. Then check if all the grids conform. Otherwise have to load them all.
        self.spacing = grid_spacing
        self.ensemble_array = None
        self.ensemble_map_dict = {}
        self.dimensions = None
        self.hist_out_dir = ""

    def relu(self, num):
        if num > 0:
            return int(num)
        else:
            return 0

    def pad_array(self, array, origin_diff, far_diff):
        """
        Takes a numpy array and pads it with the number of points specified by origin_diff
        :param array: 
        :param origin_diff: 
        :param far_diff:
        :return: numpy array
        """
        # Padding both arrays to make them the same size:
        arr = np.pad(array, ((self.relu(origin_diff[0]), self.relu(-far_diff[0])),
                           (self.relu(origin_diff[1]), self.relu(-far_diff[1])),
                           (self.relu(origin_diff[2]), self.relu(-far_diff[2]))), "constant", constant_values=0)
        return array

    def get_ensemble_array(self):
        """
        Reads in grids, converts them to 3d numpy arrays, and stacks them into 4d numpy array, which holds the information
        for the ensemble.
        :return: 
        """
        # Initialise the array
        ensemble_array = None
        # Needed for converting between Cartesian coordinates and indices.
        rec_spacing = 1.0 / self.spacing

        # Fill in ensemble array; i counts the number of grids that have been added.
        i = 0
        for p in self.paths:
            # Load in grid
            g = Grid.from_file(p)

            # Check the spacing of the grid. If different, continue and add to log.
            if g.spacing != self.spacing:
                print("Grid at {} has wrong spacing (found {}, expected {})".format(p, g.spacing, self.spacing))
                continue
            # Counter i keeps track of how many grids we have added
            i += 1

            # Get the dimensions of the grid:
            curr_dims = np.array(g.bounding_box)
            # Convert to numpy array
            arr = g.get_array()

            # Create the ensemble
            if i == 1:
                # Store the dimensions of the ensemble
                ens_dims = curr_dims
                # Put in as first element of the ensemble array
                ensemble_array = arr

            elif i == 2:
                origin_diff = (curr_dims[0] - ens_dims[0])* rec_spacing
                far_diff = ((curr_dims[1] - ens_dims[1]))* rec_spacing

                # Padding both arrays to make them the same size (and so stackable):
                arr = self.pad_array(arr, origin_diff, far_diff)
                ensemble_array =  self.pad_array(ensemble_array, -origin_diff, -far_diff)
                # Stacking 2 3D arrays creates a 4D array.
                ensemble_array = np.stack((ensemble_array, arr), axis=-1)
                # Update the ensemble dimensions
                ens_dims[0] = np.minimum(ens_dims[0], curr_dims[0])
                ens_dims[1] = np.maximum(ens_dims[1], curr_dims[1])

            else:
                origin_diff = (curr_dims[0] - ens_dims[0])* rec_spacing
                far_diff = ((curr_dims[1] - ens_dims[1]))* rec_spacing

                # Padding both arrays to make them the same size:
                arr = self.pad_array(arr, origin_diff, far_diff)
                # Ensemble array is now 4D.
                ensemble_array = np.pad(ensemble_array, ((self.relu(-origin_diff[0]), self.relu(far_diff[0])),
                                                         (self.relu(-origin_diff[1]), self.relu(far_diff[1])),
                                                         (self.relu(-origin_diff[2]), self.relu(far_diff[2])), (0, 0)),
                                        "constant", constant_values=0)
                # Np.stack stacks along a new axis, but ensemble_array is alreasy 4D, so use np.append instead.
                # When using np.append, arrays have to be the same dimension, so we expand arr with an empty 4th dimension.
                arr = np.expand_dims(arr, axis=3)
                print(arr.shape, ensemble_array.shape)
                ensemble_array = np.append(ensemble_array, arr, axis=3)
                # Update the ensemble dimensions
                ens_dims[0] = np.minimum(ens_dims[0], curr_dims[0])
                ens_dims[1] = np.maximum(ens_dims[1], curr_dims[1])

            self.dimensions = ens_dims
            self.ensemble_map_dict[i-1]= p
            self.ensemble_array = ensemble_array

    def save_gridensemble(self, path):
        """
        Pickles the GridEnsemble and saves to the location indicated in 'path'.
        :param str path:
        :return: 
        """
        pickle.dump(self, open(path), 'wb')

    def save_grid(self, array):
        """
        Given an array, outputs a grid with the dimensions of the GridEnsemble
        :param array: 3D numpy array, usually containing processed ensemble data
        :return: a :class: 'ccdc.utilities.Grid' instance
        """
        # Initialise the Grid
        grid = Grid(origin=tuple(self.dimensions[0]),
                    far_corner=tuple(self.dimensions[1]),
                    spacing=self.spacing,
                    default=0.0,
                    _grid=None)
        # Get the nonzero indices and values of the array
        nonz = array.nonzero()
        values = array[nonz]
        # Get indices per value
        as_triads = zip(*nonz)

        # Fill in the grid
        for (i, j, k), v in zip(as_triads, values):
            grid._grid.set_value(int(i), int(j), int(k), v)
        return grid

    def _plot_KS_histogram(self, dist1, dist2, tstring=""):
        """
        Plot histograms of the 2 distributions. 
        :param dist1: 1D numpy array
        :param dist2: 1D numpy array
        :param str tstring: Optional, title
        :return: 
        """
        plt.hist(dist1, bins=30, color='b', alpha=0.8, edgecolor='k')
        plt.hist(dist2, bins=30, color='r', alpha=0.5, edgecolor='k')
        plt.title(tstring)
        plt.savefig(join(self.hist_out_dir, "{}.png".format(tstring)))
        plt.close()

    def get_KS_scores(self, in_vals, threshold=2, plot=False):
        """
        Generally used for comparing ensembles (need multipleGridEnsemble?)
        in_vals has the values for one distribution as positive values, and for the other as negative.
        The function performs the Kolmogorov-Smirnov two sample test, as implemented in scipy.stats. Given two samples
        (in this case the hotspot scores from the two ensembles), the p-values gives the probability that they were sampled
        from the same distribution.(P-value close to 1- ensembles similar at that point; p-value close to zero: different.)
        :param in_vals: 1d np.array 
        :param int threshold: minimum number of values from each distribution to calculate the test on
        :param bool plot: whether to plot the distributions.
        :return: float, p_value for the KS test
        """
        positives = in_vals[np.where(in_vals > 0)]
        negatives = in_vals[np.where(in_vals < 0)]

        if len(positives) > threshold and len(negatives) > threshold:
            ks = ks_2samp(positives, -negatives)
            fstr = "KS stat: {}, pval: {}".format(ks[0], '%.3E'%Decimal(ks[1]))
            p_val = -np.log10(ks[1])
            #d = ks[0]
        else:
            fstr = "Only 1 distribution observed"
            # Hotspot scores under 10 considered not druggable/ very interesting.
            if np.all(in_vals[in_vals.nonzero()] < -10.0) or np.all(in_vals[in_vals.nonzero()] > 10.0):
                p_val = -np.log10(10 ** (-150))
                # d = 1
            else:
                p_val = 0
                # d = 0

        if plot:
            self._plot_KS_histogram(positives, -negatives, fstr)
        return p_val

    def get_2_KS_scores(self, in_vals, plot=False):
        """
        Generally used for comparing ensembles (need multipleGridEnsemble?)
        in_vals has the values for one distribution as positive values, and for the other as negative.
        The function performs the Kolmogorov-Smirnov two sample test, as implemented in scipy.stats. Given two samples
        (in this case the hotspot scores from the two ensembles), the p-values gives the probability that they were sampled
        from the same distribution.(P-value close to 1- ensembles similar at that point; p-value close to zero: different.)
        :param in_vals: 1d np.array 
        :param int threshold: minimum number of values from each distribution to calculate the test on
        :param bool plot: whether to plot the distributions.
        :return: float, p_value for the KS test
        """
        positives = in_vals[0]
        negatives = in_vals[1]

        ks = ks_2samp(positives, negatives)
        fstr = "KS stat: {}, pval: {}".format(ks[0], '%.3E'%Decimal(ks[1]))
        p_val = -np.log10(ks[1])
        # d = ks[0]

        if plot:
            self._plot_KS_histogram(positives, negatives, fstr)
        return p_val

    @staticmethod
    def plot_cluster(array):
        x, y, z = array.nonzero()
        vals = array[np.nonzero(array)]
        vals = vals.astype(np.float64)

        clust_mask = (array>0)

        unique_labels = set(vals)
        colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for k, col in zip(unique_labels, colors):
            k_mask = (array==k)
            if k == -1:
                # Black used for noise.
                col = [0, 0, 0, 1]
                noise_mask = np.where(k_mask & ~clust_mask)
                ax.scatter(noise_mask[0], noise_mask[1], noise_mask[2], c=col, s=4)
            else:
                plot_idxs = np.where(k_mask & clust_mask)
                ax.scatter(plot_idxs[0], plot_idxs[1], plot_idxs[2], c=col, s=14)
                # ax.text(plot_idxs[0][0], plot_idxs[0][1], plot_idxs[0][2], '%s' % (str(k)), size=10, zorder=1, color='k')
        plt.show()


    @staticmethod
    def DBSCAN_cluster(d_array, epsilon=4.0, mini_samples=4):
        """
        Uses ready-made clustering algorithm that infers the number of clusters
        :return: 
        """
        adj_mat = img_to_graph(d_array)
        db = DBSCAN(eps=epsilon, min_samples=mini_samples, metric="precomputed").fit(adj_mat)
        labels = db.labels_
        print(labels)

        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        print(n_clusters_)
        n_noise_ = list(labels).count(-1)
        print(n_noise_)

        a = np.zeros(d_array.shape)

        for clust, tup in zip(labels, np.ndindex(a.shape)):
            a[tup] = clust

        GridEnsemble.plot_cluster(a)
        return a

    def get_contributing_grids(self, cluster_array):
        """
        Given an array with the same first 3 dimensions of the ensemble_array, with points labelled by cluster,
        returns a dictionary of which structures contribute to which cluster.
        :param cluster_array: 3D array, labelled by cluster (eg output of self.DBSCAN_cluster())
        :return: 
        """
        clust_dic = {}
        clusters = list(set(cluster_array[cluster_array.nonzero()]))

        for c in clusters:
            cluster_mask = (cluster_array==c)
            ensemble_cluster = self.ensemble_array[cluster_mask] # should result in 2D array
            grid_indices = list(ensemble_cluster.nonzero()[1])
            clust_structs = list(set(grid_indices))
            clust_dic[c] = [(val, grid_indices.count(val)) for val in clust_structs]
        return clust_dic

    def island_find(self, arr, v, v_list=None):
        """
        Finds island within a 3d array from starting point v
        :param arr: numpy arrat
        :param v:  np.array((x, y, z)) the position of the current point
        :param v_list: the list of points within the island
        :return: list of coordinates
        """
        if not v_list:
            v_list = []
        # Points that we have visited get value of zero
        arr[tuple(v)] = 0
        # Save the point as visited
        v_list.append(tuple(v))
        # Define the neighbourhood
        radius = [2, 0, -2]
        adjacent_nodes = [np.array((x, y, z)) for x in radius for y in radius for z in radius if x or y or z]

        for node in adjacent_nodes:
            # Check that 1) point is inside array; 2) point has nonzero value
            if np.all(v + node < arr.shape) and np.all(v + node >= 0) and arr[tuple(v + node)] > 0:
                # recursively expand outwards until no points meet the criterion
                self.island_find(arr, v + node, v_list)
        # v_list now contains the coordinates of all points in the island
        return v_list

    def get_island_array(self, a):
        """
        Returns an array whose values show which cluster they belong to
        :param a: 
        :return: 
        """
        # make a copy of the array, because island_find overwrites values it's been to.
        b = np.copy(a)
        # Get the indices of all nonzero points in the array
        nonzero_a = np.nonzero(a)
        cluster_dict = {}
        counter = 1

        while len(nonzero_a[0]) > 0:
            visited_list = self.island_find(a, np.array((nonzero_a[0][0], nonzero_a[1][0], nonzero_a[2][0])), [])
            cluster_dict[counter] = visited_list
            # Update the nonzero list, because all the points in the visited cluster should have been set to zero.
            nonzero_a = np.nonzero(a)
            counter += 1

        # Returns an array whose values show which cluster they belong to.
        for key in cluster_dict.keys():
            for val in cluster_dict[key]:
                b[val] = key

        return b, cluster_dict

    def get_difference_map(self, other, tolerance=0, plot=False):
        """
        Given a second GridEnsemble with an ensemble_array of the same dimensions, performs a Kolmogorov-Smirnov test
        comparing the distributions from two 2 ensembles observed at each grid point. If tolerance is >0, values from a
        radius = tolerance from around the query point are taken into consideration. Returns a 3D numpy array, where
        the values of points correspond to the Kolmogorov-Smirnov D-value for the difference between the 2 distributions.
        
        :param other: a second GridEnsemble, whose EnsembleArray has the same first 3 dimensions as self.ensemble_array
        :param int tolerance: determines the radius around the point from which to take values into consideration.
        :return: 3D numpy array, labelled with D values for all nonzero gridpoints.
        """
        arr1 = self.ensemble_array
        arr2 = other.ensemble_array

        # Stack the two arrays along the fourth dimension. Values from the second ensemble are set to negative.
        difference_array = np.append(arr1, -arr2, axis=3)
        # Find all points that are sampled in the ensemble:
        max_diff = np.max(difference_array, axis=3) # max_diff nonzero values correspond to all points in arr1 that have been sampled.
        min_diff = np.min(difference_array, axis=3) # min_diff nonzero values correspond to all points in arr2 that have been sampled.
        diff_diff = max_diff - min_diff # nonzero indices of diff_diff indicate all points in 3D space that have hotspot values in at least 1 ensemble.

        indices = np.transpose(diff_diff.nonzero()) # get the indices per value, rather than per dimension.

        vals = difference_array[diff_diff.nonzero()] # vals is a 2D numpy array, containing the distributions of scores at each point from both ensembles.

        # idx_dict has the shape {(3D indices): Kolmogorov-Smirnov 2sample D-value}.
        idx_dict = {}
        for (a, b, c), v in zip(indices, vals):
            # Get all values within the radius specified by tolerance. Will be of shape (2*tol+1, 2*tol+1, 2*tol+1) , so flatten.
            sel = difference_array[a - tolerance: a + tolerance+1, b - tolerance: b + tolerance+1, c - tolerance: c + tolerance+1].flatten()
            # Get the Kolmogorov-Smirnov D statistic for the distributions at the sample
            d = self.get_KS_scores(sel, plot=plot)
            idx_dict[(a, b, c)] = d

        # Create an array of the modified D scores (can be used as clustering input).
        iarr = np.zeros(diff_diff.shape)
        for i, d_val in idx_dict.items():
            iarr[i] = d_val

        return iarr


if __name__ == "__main__":
    import pandas as pd
    from Ensemble import Ensemble
    from os import mkdir
    from os.path import exists, join, dirname
    import tempfile

    tempfile.tempdir = "/home/jin76872/Desktop/Mih/Data/tmp_superstar_ghecom"

    brd1_data = pd.read_csv("/home/jin76872/Desktop/Mih/Data/SIENA/EnsembleAligner/BRD1.csv")
    ref_ID = brd1_data.loc[9].squeeze()
    e = Ensemble(root_dir="", ref_ID=ref_ID)
    hot_paths = glob("/home/jin76872/Desktop/Mih/Data/SIENA/EnsembleAligner/BAZ2B/Hotspots/*/out.zip")[:100]
    hot_paths += glob("/home/jin76872/Desktop/Mih/Data/SIENA/EnsembleAligner/BRD1/Hotspots/*/out.zip")[:100]

    # trunc_paths = e.shrink_hotspots(hot_paths)
    # with open("hot_paths.txt", "w") as f:
    #     for h in trunc_paths:
    #         f.write(h + "\n")
    with open("hot_paths.txt") as f: trunc_paths=[line.strip() for line in f.readlines()]

    brd1_hots = [join(t, "acceptor.ccp4") for t in trunc_paths if "BRD1" in t]
    baz2b_hots = [join(t, "acceptor.ccp4") for t in trunc_paths if "BAZ2B" in t]

    ge1 = GridEnsemble(brd1_hots[:50])
    ge1.get_ensemble_array()
    ge1.hist_out_dir = "/home/jin76872/Desktop/Mih/Data/SIENA/saved_ensembles/acceptor_histograms_zeros_tol1"
    #ge2 = GridEnsemble(baz2b_hots)
    ge2 = GridEnsemble(baz2b_hots[:50])
    ge2.get_ensemble_array()

    # iarr = ge1.get_difference_map(ge2, tolerance=1, plot=False)
    # g_iarr = ge1.save_grid(iarr)
    # g_iarr.write(join(dirname(ge1.hist_out_dir), "BRD1_BAZ2B_50_acceptor_ranges_d_tol1.ccp4"))

    ar1 = ge1.ensemble_array
    ar2 = ge2.ensemble_array
    diffe = np.append(ar1, -ar2, axis=3)

    min_diff = np.min(diffe, axis=3)
    max_diff = np.max(diffe, axis=3)
    diff_diff = max_diff - min_diff
    vals = diffe[diff_diff.nonzero()]
    idxs = np.transpose(diff_diff.nonzero())
    pos = ar1[diff_diff.nonzero()]
    neg = ar2[diff_diff.nonzero()]

    idx_dict = {}
    for (a,b,c), p, n  in zip(idxs, pos, neg):
        sel = diffe[a-1:a+2, b-1:b+2, c-1:c+2].flatten()
        #val = ge1.get_KS_scores(v)
        val = ge1.get_2_KS_scores((p,n), plot=True)
        idx_dict[(a, b, c)] = val

    iarr = np.zeros(diff_diff.shape)
    for i, v in idx_dict.items():
        iarr[i] = v

    g_iarr = ge1.save_grid(iarr)
    g_iarr.write(join(dirname(ge1.hist_out_dir), "brd1_baz2b_acceptor_histograms_zeros_ranges_tol1.ccp4"))


    # clust_arr = GridEnsemble.DBSCAN_cluster(iarr, 5.0, 4)
    # clust_g = ge1.save_grid(clust_arr)
    # clust_g.write("BRD1_BAZ2B_diff_cluster_ranges.ccp4")
    #
    # p_vals = iarr[iarr.nonzero()]
    # nonz = iarr.nonzero()
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # p = ax.scatter(nonz[0], nonz[1], nonz[2], c=p_vals, cmap="coolwarm", alpha=0.5)
    # plt.colorbar(p)
    # plt.show()

    # diff_arr = max_diff + min_diff
    # diff_clust = GridEnsemble.DBSCAN_cluster(diff_arr)
    # clusts = list(set(diff_clust[diff_clust.nonzero()]))
    # for c in clusts:
    #     if c > 0:
    #         cvals = diffe[np.where(clust_arr == c)]
    #         cvals = cvals.flatten()
    #         ge1.get_KS_scores(cvals, plot=True)

        # nonz = diff_arr.nonzero()
    # vals = diff_arr[nonz]
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # p = ax.scatter(nonz[0], nonz[1], nonz[2], c=vals, cmap="coolwarm", alpha=0.5)
    # plt.colorbar(p)
    # plt.show()




    # max_array = np.max(ge1.ensemble_array, axis=3)
    # nonz = max_array.nonzero()
    # max_vals = max_array[nonz]
    # trunc_arr = ge1.ensemble_array[nonz]
    # freq_array = np.array([len(trunc_arr[i, :][trunc_arr[i, :] > 0]) for i in range(trunc_arr.shape[0])])
    # freq_array = freq_array * (10.0 / max(freq_array))
    # sizes = [f ** 2 for f in freq_array]
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # p = ax.scatter(nonz[0], nonz[1], nonz[2], s=sizes, c=max_vals, cmap="Greens", alpha=0.5)
    #
    # max_array = np.max(ge2.ensemble_array, axis=3)
    # nonz = max_array.nonzero()
    # max_vals = max_array[nonz]
    # trunc_arr = ge2.ensemble_array[nonz]
    # freq_array = np.array([len(trunc_arr[i, :][trunc_arr[i, :] > 0]) for i in range(trunc_arr.shape[0])])
    # freq_array = freq_array * (10.0 / max(freq_array))
    # sizes = [f ** 2 for f in freq_array]
    # # fig = plt.figure()
    # # ax = fig.add_subplot(111, projection='3d')
    # p = ax.scatter(nonz[0], nonz[1], nonz[2], s=sizes, c=max_vals, cmap="Reds", alpha=0.5)
    # fig.colorbar(p)
    # plt.show()

from __future__ import print_function, division

from ccdc.utilities import Grid
import os
from os.path import join
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy import ndimage
import pickle


class MegaGrid(object):
    """Class that handles a grid of tuples for Hotspots across multiple structures"""

    def __init__(self):
        self.stem = None
        self.prot_name = None
        self.probe = None
        self.grid_list = None
        self.tup_max_length = None
        self.results_array = None
        self.array_grid_origin = None
        self.array_grid_far_corner = None
        self.out_dir = None
        self.spacing = None

    def get_array_grid(self):
        """
        Gets the coordinates of a grid that can fit all other grids (to use as basis for self.results_array)
        :return: 
        """
        grid_list = []
        or_list = [0, 0, 0]
        far_list = [0, 0, 0]

        for root, subdirs, files in os.walk(self.stem):
            for filename in files:
                if self.probe in filename and self.prot_name in filename and 'ccp4' in filename:
                    if ('frequency' not in filename) and ('ranges' not in filename):
                        grid_list.append(join(self.stem, filename))
                        g = Grid.from_file(join(self.stem, filename))
                        _or_list = [g.bounding_box[0][j] for j in range(3)]
                        _far_list = [g.bounding_box[1][m] for m in range(3)]

                        for i in range(3):
                            or_list[i] = min(or_list[i], _or_list[i])
                            far_list[i] = max(far_list[i], _far_list[i])

        self.grid_list = grid_list
        self.spacing = g.spacing
        self.tup_max_length = len(grid_list)

        # self.array_grid = Grid(origin=(or_list[0], or_list[1], or_list[2]), far_corner=(far_list[0], far_list[1], far_list[2]),
        #               spacing=self.spacing)
        self.array_grid_origin = (or_list[0], or_list[1], or_list[2])
        self.array_grid_far_corner = (far_list[0], far_list[1], far_list[2])

    def get_results_array(self):
        """
        constructs the numpy array containing the list of tuples.
        Adjusts the coordinates based on the origins of the grids
        self.results_array = numpy array
        :return:
        """
        r_a_grid = Grid(origin=self.array_grid_origin, far_corner=self.array_grid_far_corner, spacing=self.spacing)

        nx, ny, nz = r_a_grid.nsteps
        results_array = np.zeros((nx, ny, nz), dtype=tuple)
        rec_spacing = 1 / self.spacing

        for i in range(self.tup_max_length):
            g = Grid.from_file(self.grid_list[i])
            d_coor = [g.bounding_box[0][b] - r_a_grid.bounding_box[0][b] for b in range(3)]
            print(d_coor)

            for x in range(g.nsteps[0]):
                for y in range(g.nsteps[1]):
                    for z in range(g.nsteps[2]):
                        if g.value(x, y, z) != 0:
                            r_x = x + int(d_coor[0] * rec_spacing)
                            r_y = y + int(d_coor[1] * rec_spacing)
                            r_z = z + int(d_coor[2] * rec_spacing)

                            if isinstance(results_array[r_x][r_y][r_z], tuple):
                                results_array[r_x][r_y][r_z] += (g.value(x, y, z),)
                            else:
                                results_array[r_x][r_y][r_z] = (g.value(x, y, z),)

        self.results_array = results_array

    def pickle_MegaGrid(self):
        """
        Saves MegaGrids as pickles. Always pickle the whole objects, not just the npy arrays (coordinates information lost).
        :return: 
        """
        pickle.dump(self, open(join(self.out_dir, '{}_{}_MegaGrid.p'.format(self.prot_name, self.probe)), 'wb'))

    def get_gridpoint_histograms(self):
        """
        Makes and saves histograms for each point in the results array
        """

        ind_array = np.indices(self.results_array.shape)

        def results_array_histograms(x, y, z):
            if isinstance(self.results_array[x][y][z], tuple):
                num_zeros = self.tup_max_length - len(self.results_array[x][y][z])
                if num_zeros != 0:
                    print('Num_zeros: ', num_zeros)
                hist_arr = np.array(self.results_array[x][y][z])
                # hist, bin_edges = np.histogram(hist_arr, bins=20)
                colour_dict = {"acceptor": "r", "donor": "b", "apolar": "y"}
                hist_name = self.prot_name + '_' + self.probe + '_{}_{}_{}'.format(x, y, z)

                plt.figure(1)
                plt.hist(hist_arr, bins=20, color=colour_dict[self.probe])
                plt.figtext(0.6, 0.8, ('Number of zero values:' + str(num_zeros)))
                plt.title('Score distribution at point x:{}, y:{}, z:{}'.format(x, y, z))
                plt.xlabel('Fragment hotspot score')
                plt.ylabel('Frequency')
                plt.savefig(join(self.out_dir, hist_name))
                plt.close()

        print('Generating Histograms')
        vresults_array_histograms = np.vectorize(results_array_histograms)
        vresults_array_histograms(ind_array[0], ind_array[1], ind_array[2])

    def get_gridpoint_means(self):
        """
        For each point in the 3D grid, calculates the difference in score between each point in the tuple and the mean of the tuple. 
        :return: Python list
        """
        ind_array = np.indices(self.results_array.shape)
        means = []

        def get_means(x, y, z):
            if isinstance(self.results_array[x][y][z], tuple):
                num_zeros = self.tup_max_length - len(self.results_array[x][y][z])
                if num_zeros != 0:
                    print('Hmmm')
                hist_arr = np.array(self.results_array[x][y][z])
                means.extend(list(hist_arr - np.mean(hist_arr)))

        vget_means = np.vectorize(get_means)
        vget_means(ind_array[0], ind_array[1], ind_array[2])
        return means

    def plot_gridpoint_spread(self, means):
        '''
        For each point in the 3D grid, plots the difference in score between each point in the tuple and the mean of the tuple. 
        
        '''
        mean_arr = np.array(means)
        (mu, sigma) = norm.fit(mean_arr)
        n, bins, patches = plt.hist(means, bins=40, normed=1)
        y = mlab.normpdf(bins, mu, sigma)
        a = plt.plot(bins, y, 'r--', linewidth=2)
        plt.title('Mu: {}, Sigma: {}'.format(round(mu, 2), round(sigma, 2)))
        hist_name = self.prot_name + '_{}_gridpoint_spread'.format(self.probe)
        plt.savefig(join(self.out_dir, hist_name))
        # plt.show()
        plt.close()

    def get_gridpoint_ranges(self):
        ind_array = np.indices(self.results_array.shape)
        ranges = []

        def get_ranges(x, y, z):
            if isinstance(self.results_array[x][y][z], tuple):
                num_zeros = self.tup_max_length - len(self.results_array[x][y][z])
                if num_zeros != 0:
                    print('Number of zeros: ', num_zeros)
                hist_arr = np.array(self.results_array[x][y][z])
                ranges.append(max(hist_arr) - min(hist_arr))

        vget_ranges = np.vectorize(get_ranges)
        vget_ranges(ind_array[0], ind_array[1], ind_array[2])
        return ranges

    def plot_gridpoint_ranges(self, ranges):
        '''
        Plots the range of the tuple values for each point in the 3D grid
        :return: ranges = list of the tuple ranges at each point
        '''
        ind_array = np.indices(self.results_array.shape)
        ranges = []

        def get_ranges(x, y, z):
            if isinstance(self.results_array[x][y][z], tuple):
                num_zeros = self.tup_max_length - len(self.results_array[x][y][z])
                if num_zeros != 0:
                    print('Number of zeros: ', num_zeros)
                hist_arr = np.array(self.results_array[x][y][z])
                ranges.append(max(hist_arr) - min(hist_arr))

        vget_ranges = np.vectorize(get_ranges)
        vget_ranges(ind_array[0], ind_array[1], ind_array[2])

        plt.hist(ranges, bins=40, normed=0)
        plt.title('Score ranges for {} {}'.format(self.prot_name, self.probe))
        hist_name = self.prot_name + '_{}_score_ranges'.format(self.probe)
        plt.savefig(join(self.out_dir, hist_name))
        # plt.show()
        plt.close()

    def make_frequency_grid(self):
        '''
        Makes a grid that stores how many times each point has been sampled
        :return: ccdc.utilities Grid object
        '''
        r_a_grid = Grid(origin=self.array_grid_origin, far_corner=self.array_grid_far_corner, spacing=self.spacing)

        for x in range(self.results_array.shape[0]):
            for y in range(self.results_array.shape[1]):
                for z in range(self.results_array.shape[2]):
                    if isinstance(self.results_array[x][y][z], tuple):
                        r_a_grid.set_value(x, y, z, len(self.results_array[x][y][z]))

        r_a_grid.write(join(self.out_dir, '{}_Ghecom_frequency_{}.ccp4'.format(self.prot_name, self.probe)))
        return r_a_grid

    def make_ranges_grid(self):
        '''
        Makes a grid that stores the ranges of values of each point
        :return: ccdc.utilities Grid object
        '''
        r_a_grid = Grid(origin=self.array_grid_origin, far_corner=self.array_grid_far_corner, spacing=self.spacing)

        for x in range(self.results_array.shape[0]):
            for y in range(self.results_array.shape[1]):
                for z in range(self.results_array.shape[2]):
                    if isinstance(self.results_array[x][y][z], tuple):
                        hist_arr = np.array(self.results_array[x][y][z])
                        r_a_grid.set_value(x, y, z, (max(hist_arr) - min(hist_arr)))

        r_a_grid.write(join(self.out_dir, '{}_Ghecom_ranges_{}.ccp4'.format(self.prot_name, self.probe)))
        return r_a_grid

    def make_mean_grid(self):
        '''
        Makes a grid that stores the mean of the sampld values at each point
        :return: ccdc.utilities Grid object
        '''
        r_a_grid = Grid(origin=self.array_grid_origin, far_corner=self.array_grid_far_corner, spacing=self.spacing)

        for x in range(self.results_array.shape[0]):
            for y in range(self.results_array.shape[1]):
                for z in range(self.results_array.shape[2]):
                    if isinstance(self.results_array[x][y][z], tuple):
                        t_a = np.array(self.results_array[x][y][z])
                        r_a_grid.set_value(x, y, z, np.mean(t_a))

        r_a_grid.write(join(self.out_dir, '{}_Ghecom_mean_{}.ccp4'.format(self.prot_name, self.probe)))
        return r_a_grid

    @staticmethod
    def load_MegaGrid(filename):
        """
        Loads a pickled MegaGrid
        :param filename: str, full path to pickled grid
        :return: MegaGrid object
        """
        pickle_file = open(filename, 'rb')
        newMegaGrid = pickle.load(pickle_file)

        return newMegaGrid

    def from_hotspot_maps(self, stem, out_dir, prot_name, probe_name):
        """
        Creates a MegaGrid from a number of Hotspot maps for a certain probe
        :param stem: path to directory where the grids are
        :param out_dir: path to where grids and histograms are saved
        :param prot_name: str
        :param probe_name: 'donor', 'acceptor', or 'apolar'
        :return: 
        """
        self.stem = stem
        self.out_dir = out_dir
        self.prot_name = prot_name
        self.probe = probe_name

        self.get_array_grid()
        self.get_results_array()


from __future__ import print_function, division
from hotspots.grid_extension import Grid
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.feature_extraction.image import img_to_graph
from sklearn.cluster import DBSCAN
import pickle


class GridEnsemble:
    """
    Slightly different way of compiling the Ensemble
    """

    def __init__(self, path_list, grid_spacing=0.5):
        """

        :param path_list: 
        :param grid_spacing: 
        """
        # Use the paths rather than loading in all grids at once - reduces memory usage for large ensembles.
        self.paths = path_list
        # Set the spacing before compiling the ensemble. Then check if all the grids conform. Otherwise have to load them all.
        self.spacing = grid_spacing
        self.ensemble_array = None
        self.ensemble_map_dict = {}

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
                # TODO: Make sure that we know we have not included this path for back-referencing the ensemble!
                continue
            # Counter i keeps track of how many grids we have added
            i += 1

            # Get the dimensions of the grid:
            curr_dims = np.array(g.bounding_box) * rec_spacing
            # Convert to numpy array
            arr = g.get_array()

            # Create the ensemble
            if i == 1:
                # Store the dimensions of the ensemble
                ens_dims = curr_dims
                # Put in as first element of the ensemble array
                ensemble_array = arr

            elif i == 2:
                origin_diff = (curr_dims[0] - ens_dims[0])
                far_diff = ((curr_dims[1] - ens_dims[1]))

                # Padding both arrays to make them the same size (and so stackable):
                arr = self.pad_array(arr, origin_diff, far_diff)
                ensemble_array =  self.pad_array(ensemble_array, -origin_diff, -far_diff)
                # Stacking 2 3D arrays creates a 4D array.
                ensemble_array = np.stack((ensemble_array, arr), axis=-1)
                # Update the ensemble dimensions
                ens_dims[0] = np.minimum(ens_dims[0], curr_dims[0])
                ens_dims[1] = np.maximum(ens_dims[1], curr_dims[1])

            else:
                origin_diff = (curr_dims[0] - ens_dims[0])
                far_diff = ((curr_dims[1] - ens_dims[1]))

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

            self.ensemble_map_dict[i-1]= p
            self.ensemble_array = ensemble_array

    def save_gridensemble(self, path):
        pickle.dump(self, open(path), 'wb')

if __name__ == "__main__":
    import pandas as pd
    from Ensemble import Ensemble
    from os import mkdir
    from os.path import exists

    brd1_data = pd.read_csv("/home/jin76872/Desktop/Mih/Data/SIENA/EnsembleAligner/BRD1.csv")
    ref_ID = brd1_data.loc[9].squeeze()
    e = Ensemble(root_dir="", ref_ID=ref_ID)
    hot_paths = glob("/home/jin76872/Desktop/Mih/Data/SIENA/EnsembleAligner/BAZ2B/Hotspots/*/out.zip")[:10]
    hot_paths += glob("/home/jin76872/Desktop/Mih/Data/SIENA/EnsembleAligner/BRD1/Hotspots/*/out.zip")[:10]

    trunc_paths = e.shrink_hotspots(hot_paths)
    with open("hot_paths.txt", "w") as f:
        for h in hot_paths:
            f.write(h)

    brd1_hots = [t for t in trunc_paths if "BRD1" in t]
    baz2b_hots = [t for t in trunc_paths if "BAZ2B" in t]

    ge1 = GridEnsemble(brd1_hots)
    ge2 = GridEnsemble(baz2b_hots)

    diffe = np.append(ge1.ensemble_array, -ge2.ensemble_array)

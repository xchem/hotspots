from __future__ import print_function, division
from glob import glob
from ccdc.protein import Protein
from hotspots.grid_extension import Grid, _GridEnsemble
from scipy.ndimage import morphology
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull


class CADCalculator(object):
    """
    Class to calculare the CAD matrix for a single input pdb.
    TODO: Make this work on selections (ie binding site).
    """
    def __init__(self, protein_path, grid_spacing=0.5, chains=None):
        """
        :param str protein_path: path to target protein to be CADed
        :param float grid_spacing: how far away grid points are (in Angstroms). Finer spacing gives better resolution at the cost of speed. 
        """
        self.prot_path = protein_path
        self.g_spacing = grid_spacing
        self.chains = chains
        if not chains:
            self.chains = ["A"]

        # variables that may be changed later
        self.protein = None
        self.prot_sequence = None
        self.aa_grid = None
        self.CAD_matrix = None

    def get_protein(self):
        """
        Gets the relevant bits from the pdb.
        :return: 
        """
        prot = Protein.from_file(self.prot_path)
        for chain in prot.chains:
            print(chain.identifier)
            if chain.identifier not in self.chains:
                prot.remove_chain(chain.identifier)
        prot.remove_all_waters()
        for ligand in prot.ligands:
            prot.remove_ligand(ligand.identifier)
        prot.remove_all_metals()
        prot.remove_hydrogens()

        self.protein = prot

    def get_sequence(self):
        """
        Get the sequence for the selected chains of the target protein.
        :return: 
        """
        if not self.protein:
            self.get_protein()
        self.sequence = [r.identifier for r in self.protein.residues]

    @staticmethod
    def set_uniform_values(grd, fill=1):
        nx, ny, nz = grd.nsteps
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                        if grd.value(i, j, k) != 0.0:
                            grd.set_value(i, j, k, fill)

        return grd

    def get_residue_vdw(self, residue, res_num, padding=3.0):
        """
        generates a mask for the residue where points within VdW distance of residue heavy atoms are 
        :param residue: `ccdc.prottein.Residue`
        :param res_num: The value which to assign to voxels in the masked area
        :param padding: float, padding around minimal coordinates in Angstroms
        :return: `hotspots.grid_extension.Grid`
        """
        coords = np.array([a.coordinates for a in residue.atoms])
        min_coords = (np.min(coords[:, 0]), np.min(coords[:, 1]), np.min(coords[:, 2]))
        max_coords = (np.max(coords[:, 0]), np.max(coords[:, 1]), np.max(coords[:, 2]))

        # Put some padding around the minimal and maximum values:
        g_origin = tuple(x - padding for x in min_coords)
        g_far_corner = tuple(y + padding for y in max_coords)

        layer = Grid(origin=g_origin, far_corner=g_far_corner, spacing=self.g_spacing, default=0.0, _grid=None)

        for a in residue.atoms:
            layer.set_sphere(point=a.coordinates,
                         radius=a.vdw_radius,
                         value=1,
                         scaling='None')
        layer = self.set_uniform_values(layer, res_num)
        print("Size of layer: {}".format(layer.count_grid()))

        return layer

    def grid_from_protein(self):
        """
        Puts the protein inside a CCDC Grid object
        :return: a :class: hotspots.grid_extension.Grid instance
        """
        if not self.protein:
            self.get_protein()

        coords = np.array([a.coordinates for res in self.protein.residues for a in res.atoms])
        min_coords = (np.min(coords[:, 0]), np.min(coords[:, 1]), np.min(coords[:, 2]))
        max_coords = (np.max(coords[:, 0]), np.max(coords[:, 1]), np.max(coords[:, 2]))

        # Put some padding around the minimal and maximum values:
        g_origin = tuple(x-3.0 for x in min_coords)
        g_far_corner = tuple(y + 3.0 for y in max_coords)

        prot_grid = Grid(origin=g_origin, far_corner=g_far_corner, spacing=self.g_spacing, default=0.0, _grid=None)
        return prot_grid

    def get_aa_grid(self):
        """
        Returns a grid where points within the VdW radius of e residue are labelled with the residue's index+1
        :return: 
        """
        #pg = self.grid_from_protein()
        if not self.protein:
            self.get_protein()

        res_names = [res.identifier for res in self.protein.residues]
        res_nums = [idx + 1 for (idx, val) in enumerate(res_names)]

        res_gs = [self.get_residue_vdw(r, res_nums[i]) for i, r in enumerate(self.protein.residues)]
        print(len(set(round(g.spacing, 3) for g in res_gs)))
        #res_gs = Grid.common_grid(res_gs, padding=0)
        ge = _GridEnsemble()
        self.aa_grid = ge.from_grid_list(res_gs, os.getcwd(), "NUDT5114", "acceptor")
        self.aa_grid.write(self.prot_path.replace(".pdb", "_vdw_grid.ccp4"))


    def get_contact_areas(self):
        """
        
        :return: 
        """
        if not self.aa_grid:
            self.get_aa_grid()

        vdw_array = self.aa_grid.get_array()
        res_names = [res.identifier for res in self.protein.residues]
        res_nums = [idx + 1 for (idx, val) in enumerate(res_names)]

        #Not contact dict -> directly fill in the matrix. See how slow it is later.
        cad_matrix = np.zeros((len(res_nums), len(res_nums)))

        struct_element = morphology.generate_binary_structure(rank=3, connectivity=2)
        small_struct_element = morphology.generate_binary_structure(3, 2)
        large_struct_element = morphology.iterate_structure(structure=struct_element, iterations=2)
        print(large_struct_element)

        for r in res_nums:
            loc_arr = np.zeros(vdw_array.shape)
            loc_arr[np.where(vdw_array==r)]=1
            dil = morphology.binary_dilation(loc_arr, large_struct_element, iterations=2)
            contacts = dil*vdw_array
            eroded_contacts = morphology.binary_erosion(dil, small_struct_element)
            eroded_contacts = np.abs(eroded_contacts-1)
            contacts = contacts*eroded_contacts
            near_res = list(set(contacts[contacts.nonzero()]))

            print(r, near_res)
            for n in near_res:
                if int(n) != int(r):
                    cad_matrix[int(r-1), int(n-1)] = len(contacts[np.where(contacts==n)])
                    # idx = np.where(contacts==n)
                    # if np.array(zip(*idx)).shape[0] >=4:
                    #     hull = ConvexHull(np.array(zip(*idx)))
                    #     cad_matrix[int(r - 1), int(n - 1)] = hull.area
                    # else:
                    #     cad_matrix[int(r - 1), int(n - 1)] ==1.0
            if r==1 or r==50:
                g = CADCalculator.from_array(contacts, self.aa_grid.bounding_box[0], self.aa_grid.bounding_box[1])
                g.write("contacts_g_{}_acceptor_ranges.ccp4".format(self.protein.residues[r-1]))

        self.CAD_matrix = cad_matrix
        fname = self.prot_path.replace(".pdb", "_CAD.txt")
        np.savetxt(fname, cad_matrix)

    @staticmethod
    def from_array(array, g_origin, g_far_corner):
        """
        creates a grid from array
        :param fname: path to pickled numpy array
        :return: `hotspots.grid_extension.Grid`
        """
        grid = Grid(origin=g_origin,
                    far_corner=g_far_corner,
                    spacing=0.5,
                    default=0.0,
                    _grid=None)

        indices = np.nonzero(array)
        values = array[indices]
        as_triads = zip(*indices)

        for (i, j, k), v in zip(as_triads, values):
            grid._grid.set_value(int(i), int(j), int(k), v)

        return grid

    @staticmethod
    def load_cad(filepath):
        """
        Creates a CADCalculator from a precalculated matrix
        :param filepath: 
        :return: 
        """
        cad_mat = np.loadtxt(filepath)
        prot_name = filepath.replace("_CAD.txt", ".pdb")
        new_CC = CADCalculator(prot_name)
        new_CC.CAD_matrix = cad_mat
        return new_CC

    @staticmethod
    def get_cad_score(ref_cad, model_cads):
        """
        Calculates CAD-score as descibed in Abagyan and Totrov, 1997, Eq.1
        
        :param ref_cad: CADCalculator object for the reference model (with cad_mat)
        :param model_cads: CADCalculator objects for the query structures
        :return: dictionary of model name: CAD_score
        """
        # Values taken from the paper:
        c = 0.9
        B_std = 20.0

        # Assume B-factors == 1 for all residues (same as paper).
        b_factor_i = 1.0
        b_factor_j = 1.0
        wi = np.exp(- b_factor_i / B_std)
        wj = np.exp(- b_factor_j / B_std)

        ref_mat = ref_cad.CAD_matrix

        score_dict = {}
        for m_cad in model_cads:
            diff_mat = abs(ref_mat - m_cad.CAD_matrix)
            sum_mat = ref_mat + m_cad.CAD_matrix
            Aworst = c * np.sum(wi * wj * 0.5 * sum_mat)
            score = 100.0 * np.sum(wi * wj * diff_mat) / Aworst
            score_dict.update({os.path.basename(m_cad.prot_path): score})

        return score_dict

if __name__ == "__main__":
    # p_path = glob(os.path.join(os.getcwd(), "NUDT5", "*bound*.txt"))
    ref_struct = glob(os.path.join(os.getcwd(), "*ground*.pdb"))[0]
    ref_cc = CADCalculator(protein_path=ref_struct, chains=["A", "B"])
    #ref_cc = CADCalculator(protein_path=ref_struct)#
    ref_cc.get_contact_areas()
    #ref_struct = p_path[-3]
    #ref_cc = CADCalculator.load_cad(ref_struct)

    # p_path = glob(os.path.join(os.getcwd(), "NUDT5", "*bound*.txt"))
    # cc_list = [CADCalculator.load_cad(p) for p in p_path]
    # scores = CADCalculator.get_cad_score(ref_cc, cc_list)

    p_path = glob(os.path.join(os.getcwd(), "NUDT5", "*bound*.pdb"))[0:10]
    cad_mats = []
    for pp in p_path:
        cc = CADCalculator(protein_path=pp, chains=["A", "B"])
        cc.get_contact_areas()
        cad_mats.append(cc)
        # del(cc)
    scores = CADCalculator.get_cad_score(ref_cc, cad_mats)

    # for m in range(len(cad_mats)):
    #     plt.imshow(cad_mats[m], cmap="hot", interpolation="nearest")
    #     plt.savefig(os.path.join(os.getcwd(), "cad_{}.png").format(str(m)))
    #     plt.close()


    # ref_cc.get_contact_areas()
    # plt.imshow(ref_cc.CAD_matrix, cmap="hot", interpolation="nearest")
    # plt.savefig(os.path.join(os.getcwd(), "cad_ground_state.png"))
    # plt.close()


    # cad_scores = []
    # for m in cad_mats:
    #     wi = np.exp(-1/20)
    #     wj = np.exp(-1/20)
    #     C = 0.9
    #     diff_mat = ref_cc.CAD_matrix - m
    #     sum_mat = ref_cc.CAD_matrix +m
    #     Aworst = C*np.sum(wi*wj*0.5*abs(sum_mat))
    #     #vals = diff_mat[diff_mat.nonzero()]
    #     #score = (sum(vals), sum(abs(vals)))
    #     score = 100.0*np.sum(wi*wj*abs(diff_mat))/Aworst
    #     cad_scores.append(score)




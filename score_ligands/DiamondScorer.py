from ccdc.protein import Protein
import os
from os.path import exists
from glob import glob
from hotspots import hs_io, grid_extension
from hotspots.hs_utilities import Helper
from ccdc.io import MoleculeWriter, MoleculeReader
import numpy as np
import datetime
import shutil


class DiamondScorer(object):
    """
    Class to score fragments in each hotspot. Works with individual hotspots results.
    """
    def __init__(self, stem, prot_name, x_id, keyword=None, out_dir=None, chains=None, data_dir=None):
        """
        
        :param str stem: directory holding hotspot results for different structures of the target. 
        :param str prot_name: name of the protein (eg NUDT5)
        :param str x_id: the crystal id (eg "x144")
        :param str keyword: Any additional identifier needed: eg "bound" or "ground" state of PanDDA model.
        :param str out_dir: ouput directory for the scored ligands. Defaults to data_dir.
        :param list chains: In cases where hotspot maps have been calculated separately for different chains of the protein, 
                        specifies which result to use. Defaults to ["A"].
        :param str data_dir: dir where our target hotspot result lives
        """
        self.keyword = keyword
        self.protein_name = prot_name
        self.xstal_id = x_id
        self.out_dir = out_dir
        self.chains = chains
        if not self.chains:
            self.chains = ["A"]
        self.data_dir = data_dir
        if not self.data_dir:
            self.get_data_dir(stem)
        if not self.out_dir:
            self.out_dir= self.data_dir


    def get_data_dir(self, stem):
        """
        
        :return: 
        """
        #res_dir = os.path.join(self.stem, "hotspot_results_try_1")
        #assert exists(res_dir), "Cannot find results directory for the target"
        str_chains = "".join(self.chains)
        if self.keyword:
            out = os.path.join(stem, "{}_{}_{}_{}_results".format(self.protein_name, self.xstal_id, str_chains, self.keyword))

        else:
            out = os.path.join(stem, "{}_{}_{}_results".format(self.protein_name, self.xstal_id, str_chains))
            if not exists(out):
                print("Can't find xstal-specific results dir")

        self.data_dir = out

    def get_hotspot(self):
        """
        
        :return: 
        """
        hs_reader = hs_io.HotspotReader(os.path.join(self.data_dir, "out.zip"))
        hs_result = hs_reader.read()
        return hs_result

    def get_ligands(self):
        """
        
        :return: 
        """
        lig_reader = MoleculeReader(os.path.join(self.data_dir, "ligands.mol2"))

        return lig_reader

    def read_protein(self, hotspot_result=None):
        """
        Reads in the protein used to calculate the hotspot maps in data_dir.
        :return: 
        """
        if not hotspot_result:
            hr = self.get_hotspot()
        else:
            hr = hotspot_result

        return hr.protein

    def get_ligand_chain(self, ligand, hotspot_result=None):
        """
        For each ligand, determine which protein chain (or chains) it is closest to.
        :param ccdc.molecule.Molecule ligand: 
        :return: 
        """
        prot = self.read_protein(hotspot_result)
        bs = Protein.BindingSiteFromMolecule(prot, ligand, 5)
        ligand_chains = list(set([r.identifier.split(":")[0] for r in bs.residues]))
        return ligand_chains

    def get_ligand_features(self, ligand, hotspot_result=None):
        """
        Checks if fragment falls in favourable area of the hotspot map. Assigns that area to the three 
        closest amino acids.
        :param ligand: 
        :param hotspot_result: 
        :return: 
        """
        if not hotspot_result:
            hotspot_result = self.get_hotspot()

        atoms = ligand.heavy_atoms
        feats = []

        for a in atoms:
            if a.is_donor and a.is_acceptor:
                tar_grids = {"donor": hotspot_result.super_grids["donor"],
                             "acceptor": hotspot_result.super_grids["acceptor"]}
            elif a.is_donor:
                tar_grids = {"donor": hotspot_result.super_grids["donor"]}

            elif a.is_acceptor:
                tar_grids = {"acceptor": hotspot_result.super_grids["acceptor"]}
            else:
                continue

            for key, g in tar_grids.items():
                #tar_isls = [isl for isl in g.islands(1) if isl.contains_point(a.coordinates, tolerance=1.0)]
                tar_isls = [isl for isl in g.islands(1) if isl.value_at_coordinate(a.coordinates, tolerance=2)[0] > 0.0]
                if len(tar_isls) != 0:
                    feat = LigandFeature(scored_atom=a, scored_atom_type=key, feature_grids=tar_isls)
                    feat.get_proxy_residues(hotspot_result.protein)
                    feats.append(feat)

        return feats

    @staticmethod
    def scaled_value_at_coordinate(grid, coordinates, tolerance=2):
        """
        Uses Grid.value() rather than Grid.value_at_point(). Chris Radoux reported speed issues.
        :param coordinates:
        :param tolerance:
        :return:
        """
        i, j, k = grid.point_to_indices(coordinates)
        nx, ny, nz = grid.nsteps
        scores = {}
        # Calculating the maximum radius (in Angstroms) where points will be scored
        radius = 3.0**(0.5)*grid.spacing*tolerance + 0.5*3.00**(0.5)*grid.spacing
        recip_radius = radius**(-1)

        for di in range(-tolerance, +tolerance + 1):
            for dj in range(-tolerance, +tolerance + 1):
                for dk in range(-tolerance, +tolerance + 1):
                    if 0 <= (i + di) < nx and 0 <= (j + dj) < ny and 0 <= (k + dk) < nz:
                        scores.update({grid.value(i + di, j + dj, k + dk): (i + di, j + dj, k + dk)})

        if len(scores) > 0:
            scaled_scores={}
            for score, gpoint in scores.items():
                if score < 0.1:
                    sc_score = 0
                    point = (0, 0, 0)

                else:
                    a, b, c = gpoint
                    point = grid.indices_to_point(a, b, c)
                    distance = Helper.get_distance(coordinates, point)
                    sc_score = (radius-distance)*recip_radius*score
                scaled_scores.update({sc_score: point})

            max_scaled_score = sorted(scaled_scores.keys(), reverse=True)[0]
            max_pos = scaled_scores[max_scaled_score]


        else:
            max_scaled_score = 0
            max_pos = (0, 0, 0)

        return max_scaled_score, max_pos

    def get_scaled_score(self, ligand, tolerance=2, hotspot=None):
        """
        
        :param hotspot: 
        :param tolerance: 
        :return: 
        """
        if not hotspot:
            hotspot = self.get_hotspot()
        lig_scores = {}
        for a in ligand.heavy_atoms:
            if a.is_donor and a.is_acceptor:
                score_dict = {}
                acc_score, acc_pos = self.scaled_value_at_coordinate(hotspot.super_grids["acceptor"], a.coordinates, tolerance)
                don_score, don_pos = self.scaled_value_at_coordinate(hotspot.super_grids["donor"], a.coordinates, tolerance)
                score_dict.update({acc_score: acc_pos,
                                   don_score: don_pos})
                max_score = sorted(score_dict.keys(), reverse=True)[0]
                max_pos = score_dict[max_score]

            elif a.is_acceptor:
                max_score, max_pos = self.scaled_value_at_coordinate(hotspot.super_grids["acceptor"], a.coordinates, tolerance)

            elif a.is_donor:
                max_score, max_pos = self.scaled_value_at_coordinate(hotspot.super_grids["donor"], a.coordinates, tolerance)

            else:
                max_score, max_pos = self.scaled_value_at_coordinate(hotspot.super_grids["apolar"], a.coordinates,
                                                                     tolerance)
                max_score=0

            lig_scores.update({a:(max_score, max_pos)})
            a.partial_charge = max_score

        return ligand, lig_scores


    def scaled_score_ligands(self, tolerance):
        """
        Applies linear scaling to scores assigned to atom, depending on distance between atom and the scored point.
        :param int tolerance: How many gridpoints away is it acceptable for an atom to be from the nearest point of its corresponding map.
        :return: 
        """
        hs = self.get_hotspot()
        all_ligs = self.get_ligands()

        ligs = {}
        # score only ligands from chains where we have hotspots info
        for al in all_ligs:
            lchains = self.get_ligand_chain(ligand=al, hotspot_result=hs)
            if set(lchains) == set(self.chains):
                ligs[al] = lchains

        if len(ligs) == 0:
            print("No ligands to score")
            return

        if not self.out_dir:
            self.out_dir = self.data_dir

        scored_ligs = []

        for lig, lig_chains in ligs.items():
            scored_lig = self.get_scaled_score(lig, tolerance, hs)[0]
            joined_lig_chains = "".join(lig_chains)
            ligand_score = np.mean([a.partial_charge for a in scored_lig.heavy_atoms])
            scored_lig.identifier += "_{}_{}_{}_{}".format(self.protein_name, self.xstal_id, joined_lig_chains, round(ligand_score, 2))
            scored_ligs.append(scored_lig)


        with MoleculeWriter(os.path.join(self.out_dir, "scaled_scored_ligands.mol2")) as writer:
            for ligand in scored_ligs:
                writer.write(ligand)


    def score_ligands(self):
        """
        Scores ligand using the hotspots API scoring function
        :return: 
        """
        hs = self.get_hotspot()
        all_ligs = self.get_ligands()

        ligs = {}
        # score only ligands from chains where we have hotspots info
        for al in all_ligs:
            lchains = self.get_ligand_chain(ligand=al, hotspot_result=hs)
            if set(lchains) == set(self.chains):
                ligs[al] = lchains

        if len(ligs) == 0:
            print("No ligands to score")
            return

        if not self.out_dir:
            self.out_dir = self.data_dir

        scored_ligs = []

        for lig, lig_chains in ligs.items():
            scored_lig = hs.score(lig)
            # change the identifier of the scored_lig to include the protein chain it came from
            joined_lig_chains = "".join(lig_chains)
            ligand_score = np.mean([a.partial_charge for a in scored_lig.heavy_atoms])
            scored_lig.identifier += "_{}_{}_{}_{}".format(self.protein_name, self.xstal_id, joined_lig_chains, round(ligand_score, 2))
            scored_ligs.append(scored_lig)

        with MoleculeWriter(os.path.join(self.out_dir, "scored_ligands.mol2")) as writer:
            for ligand in scored_ligs:
                writer.write(ligand)

        #return scored_ligs

    def log_DiamondScorer(self):
        """
        Produces a log of the scored ligands
        :return: 
        """
        if not self.out_dir:
            self.out_dir = self.data_dir

        now = datetime.datetime.now()

        with open(os.path.join(self.out_dir, "ligand_scoring.log"), "a") as f:
            f.write("DiamondScorer log calculated at {}-{}-{} {}:{} \n".format(now.year, now.month, now.day, now.hour,
                                                                          now.minute))
            for key, val in self.__dict__.items():
                f.write(str(key) + ": " + str(val) + "\n")


class LigandFeature(object):
    """
    Stores information for areas of the hotspot maps sampled by the fragments
    """
    def __init__(self, scored_atom, scored_atom_type, feature_grids):
        self.atom = scored_atom
        self.feature_type = scored_atom_type
        self.grids = feature_grids
        self.atom_score = [i.value_at_coordinate(self.atom.coordinates, tolerance=2) for i in self.grids]
        self.feature_size = [g.count_grid() for g in self.grids]
        self.centroids = [g.centroid() for g in self.grids]
        self.distances = [Helper.get_distance(self.atom.coordinates, self.atom_score[i][1]) for i in range(len(self.atom_score))]
        self.centroid_proxy_residues = None
        self.scored_atom_proxy_residues = None


    @staticmethod
    def get_close_residues(coordinates, protein, max_distance, atom_type):
        """
        Picks out the closest 3 residues that could contribute to a feature
        :param coordinates: 
        :param protein: 
        :param float max_distance: what radius to look for residues
        :return: 
        """
        if atom_type == "donor":
            atoms = [a for a in protein.atoms if a.is_acceptor]
        elif atom_type == "acceptor":
            atoms = [a for a in protein.atoms if a.is_donor]
        else:
            print("Unsupported interaction type in LigandFeature.get_close_residues()")
            return

        near_atoms = {}
        for atm in atoms:
            dist = Helper.get_distance(atm.coordinates, coordinates)
            if dist < max_distance:
                if dist in near_atoms.keys():
                    near_atoms[dist].append(atm)
                else:
                    near_atoms.update({dist: [atm]})
            else:
                continue
        if len(near_atoms.keys()) == 0:

            return None

        else:
            #print(near_atoms)
            if len(near_atoms.keys()) >3:
                closest = sorted(near_atoms.keys())[0:3]
            else:
                closest = near_atoms.keys()

            select = [near_atoms[c] for c in closest]
            #print(select)

            return [(m.label, m.residue_label) for s in select for m in s]
            #return [m for s in select for m in s]

    def get_proxy_residues(self, hs_protein):
        """
        
        :param protein: a :class: ccdc.protein.Protein instance
        :return: 
        """

        self.centroid_proxy_residues = [self.get_close_residues(coordinates=p,
                                                                   protein=hs_protein,
                                                                   max_distance=5.0,
                                                                   atom_type=self.feature_type) for p in self.centroids]

        self.scored_atom_proxy_residues = self.get_close_residues(coordinates=self.atom.coordinates,
                                                                  protein=hs_protein,
                                                                  max_distance=5.0,
                                                                  atom_type=self.feature_type)
        print(self.scored_atom_proxy_residues)



    def get_line_of_sight(self):
        """
        
        :return: 
        """
        for a in self.scored_atom_proxy_residues:
            print(self.atom, a, a.is_in_line_of_sight(self.atom))












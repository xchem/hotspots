from __future__ import print_function, division
from ccdc.protein import Protein
from ccdc.io import MoleculeWriter, MoleculeReader
import pandas as pd
from os.path import join, basename, dirname, exists
from os import mkdir
from glob import glob
from hotspots import data
from hotspots.hs_io import HotspotReader, HotspotWriter
from hotspots.grid_extension import Grid
import numpy as np
from GridEnsemble import GridEnsemble
import luigi
from run_hotspots import ParalleliselRunner


class Ensemble:
    """
    Parent class with methods used by both EnsembleResult and MultipleEnsemble
    """
    def __init__(self, root_dir=None, ref_ID=None):
        # The pipeline root directory
        self.root_dir = root_dir
        # Pandas series, carries information on the reference structure for the ensemble
        self.reference_ID = ref_ID
        # The extracted binding site of the reference structure (ccdc.Protein.BindingSite)
        self.reference_binding_site = None
        # Although the Ensemble class doesn't do any alignments, the child classes both use the same settings.
        self.alignment_settings = Protein.ChainSuperposition.Settings()
        # Set the default for this from a config?
        self.alignment_settings.superposition_atoms = "BACKBONE"
        # TODO: controlled by config file - linked to the hotspot settings
        self.hotspot_probe_types = ["donor", "acceptor", "apolar", "positive", "negative"]

    @staticmethod
    def get_binding_site(row):
        """
        Given a row from an EnsembleResult.protein_data dataframe, extract the binding site of interest
        :param row: 
        :return: a :class: ccdc.protein.BindingSiteFromListOfResidues instance
        """
        # Load up protein
        print("Getting binding site for {}".format(basename(row["Filename"])))
        p = Protein.from_file(row["Filename"])

        # Remove all chains we don't care about
        print(row["Chains"])
        for chain in p.chains:
            if chain.identifier not in row["Chains"]:
                p.remove_chain(chain.identifier)

        # Get correct bonding for the ligand
        p.detect_ligand_bonds()

        # Eliminate ligands from other chains and binding sites.
        for lig in p.ligands:
            # Get the residues around the ligand (5 A radius).
            lbs = Protein.BindingSiteFromMolecule(p, lig, 5)
            try:
                ligand_chains = list(set([r.identifier.split(":")[0] for r in lbs.residues]))

            # CCDC API can fail here if the ligand is an aminoacid (eg in bromodomains), as it then looks for the amino acid in the protein sequence.
            except IndexError:
                continue
            # If there is no overlap between the ligand chain and the target chains, remove the ligand.
            if all(chain not in row["Chains"] for chain in ligand_chains):
                p.remove_ligand(lig.identifier)

        # Select the binding site residues
        resi = [r for r in p.residues if r.identifier in row["Binding site"]]

        # Create the binding site object
        bs = Protein.BindingSiteFromListOfResidues(p, resi)
        # Set the identifier of the protein for easy recall later
        bs.protein.identifier = "{}_{}".format(row["Ensemble ID"], row["ID"])
        print("Chains", [ch.identifier for ch in bs.protein.chains])
        return bs

    @staticmethod
    def shrink_to_binding_site(in_grid, new_origin, new_far_corner):
        """
        Given an input grid, will reduce it to the area defined by the new origin and far corner
        :param in_grid: a ccdc.utilities.Grid
        :param new_origin: numpy array((x, y, z))
        :param new_far_corner: numpy array((x, y, z))
        :return: ccdc.utilities.Grid
        """
        # Check that the new coordinates fall within the grid:
        ori = np.array(in_grid.bounding_box[0])
        far_c = np.array(in_grid.bounding_box[1])

        if (new_origin > ori).all() and (new_far_corner < far_c).all():
            ori_idxs = in_grid.point_to_indices(tuple(new_origin))
            far_idxs = in_grid.point_to_indices(tuple(new_far_corner))
            region = ori_idxs + far_idxs
            # Get only the sub-grid defined by the 6 indices in region
            new_g = in_grid.sub_grid(region)
            return new_g

        else:
            print("Selected area larger than grid; try reducing the padding in shrink_hotspots()")
            # TODO: Log as error

    def shrink_hotspots(self, hotspot_paths, padding=2.0):
        """
        Takes in the calculated hotspots on the aligned ensemble. Crops and saves only the area around the reference binding site.
        Results are stored in the same parent directory as the fullsized hotspots, in dir called "binding_site_maps"
        :param list hotspot_paths: Paths to the hotspot results we would like to shrink.
        :param float padding: How many angstroms away from furthest binding site atom to look.
        :return: list of the paths for all shrunk hotspots in the ensemble.
        """
        # Get the area to truncate around the binding site:
        print("Shrinking hotspots for ensemble...")

        if not self.reference_binding_site:
            self.reference_binding_site = self.get_binding_site(self.reference_ID)

        # Find the maximum and minimum coordinates of the reference binding site
        dims = np.array([a.coordinates for r in self.reference_binding_site.residues for a in r.atoms])
        min_coords = np.array([np.min(dims[:, 0]), np.min(dims[:, 1]), np.min(dims[:, 2])])
        max_coords = np.array([np.max(dims[:, 0]), np.max(dims[:, 1]), np.max(dims[:, 2])])

        # Add some padding in both directions:
        min_coords -= padding
        max_coords += padding

        h_out_dir_list = []

        for p in hotspot_paths:
            # Read in hotspot result
            h_result = HotspotReader(p).read()
            # Shrink the grids for each probe type
            for probe, g in h_result.super_grids.items():
                h_result.super_grids[probe] = self.shrink_to_binding_site(g, min_coords, max_coords)

            res_path = dirname(p)
            # Save shrunk hotspot, assuming the directory it was previously in was named sensibly.
            h_out_dir = join(res_path, "binding_site_maps")
            h_out_dir_list.append(join(h_out_dir, "out"))
            with HotspotWriter(h_out_dir, visualisation="pymol", grid_extension=".ccp4", zip_results=False) as writer:
                writer.write(h_result)

        return h_out_dir_list

class EnsembleResult(Ensemble):
    """
    Class for storing ensembles of structures of a single target. Also deals with bound ligands.
    """

    def __init__(self, root_dir, ref_id, df, aligned=False):
        Ensemble.__init__(self, root_dir, ref_id)
        self.ensemble_ID = self.reference_ID["Ensemble ID"]  # Usually the target name
        self.protein_data = df  # pd Dataframe, output of SIENAReader or XChemReader
        self.ligand_data = None  # pd Dataframe, based on the ligands extracted from structures in protein_data
        self.aligned = aligned  # whether the structures in the ensemble are aligned or not.

    def _save_protein(self, bs):
        """
        Saves the preprocessed, aligned proteins, each in its own folder.
        :param bs: 
        :return: str, the path to where the protein is saved
        """
        print("Saving protein", bs.protein.identifier, [ch.identifier for ch in bs.protein.chains])
        p = bs.protein
        splitid = p.identifier.split("_")
        assert len(splitid) == 2, "Identifier not in the right format"
        fname = p.identifier + ".pdb"

        prot_dir = join(self.root_dir, p.identifier)

        if not exists(prot_dir):
            mkdir(prot_dir)

        p_fname = join(prot_dir, fname)
        with MoleculeWriter(p_fname) as protein_writer:
            protein_writer.write(p)

        return p_fname

    def align_ensemble(self, bs_list):
        """
        Aligns all elements of the ensemble to the reference structure
        :param list bs_list: python list of :class: 'ccdc.Protein.BindingSiteFromListOfResidues' instances
        :return: list of aligned :class: 'ccdc.Protein.BindingSiteFromListOfResidues instances'
        """
        print("Aligning ensemble proteins for {}".format(self.ensemble_ID))
        # Keep track of the RMSDs between the reference and the ensemble members
        rmsd_list = []
        # Copy the binding site to avoid changing it during superposition - investigate if actually needed.
        e_ref_bs = self.reference_binding_site

        for bs in bs_list:
            chain_superposition = Protein.ChainSuperposition(self.alignment_settings)
            (rmsd, transformation) = chain_superposition.superpose(e_ref_bs.protein.chains[0],
                                                                   bs.protein.chains[0],
                                                                   e_ref_bs)
            rmsd_list.append(rmsd)

        # Add RMSD column to the protein dataframe
        self.protein_data["RMSD"] = rmsd_list
        self.aligned= True
        return bs_list

    def process_ensemble_proteins(self):
        """
        If proteins need further alignment - does that. Saves only relevant chains and ligands within radius of the target chain.
        :return: 
        """
        # Load up the binding sites of the ensemble:
        bs_list = [self.get_binding_site(r) for idx, r in self.protein_data.iterrows()]

        if not self.aligned:
            bs_list = self.align_ensemble(bs_list)
            self.aligned = True

        paths_list = [self._save_protein(bs) for  bs in bs_list]
        self.protein_data.to_csv(join(self.root_dir, "{}.csv".format(self.ensemble_ID)))
        return paths_list

    def run_hotspots(self, paths_list, nrotations=100000, charged=False):
        """
        :param int nrotations: how many probe rotations to do in the hotspots calculation. Supply in config file.
        :param bool charged: whether to use charged probes in the hotspots calculation. Supply in config file.
        :return: return the paths to the hotspots.
        """

        # Not the most intelligent use of Luigi, but will calculate a lot of hotspots very reliably.
        luigi.build([ParalleliselRunner(paths_list, nrotations, charged)], local_scheduler=True, workers=30)
        # TODO: see what happens in case luigi fails on some of these
        hotspot_paths_list = [join(dirname(p), "fullsize_hotspots", "out.zip") for p in paths_list]

        return hotspot_paths_list

    def make_ensemble_maps(self, hotspot_rotations, hotspot_charged_probes):
        """
        Make an ensemble map for each available probe type and save it.
        :return: 
        """
        hots_list = self.run_hotspots(hotspot_rotations, hotspot_charged_probes)

        # This returns the directories where the shrunk hotspots are stored, but not paths to individual maps.
        shrunk_hots = self.shrink_hotspots(hots_list)

        # TODO Make sure that GridEnsemble can handle not being supplied paths in a more graceful way.
        for probe in self.hotspot_probe_types:
            ge_paths = glob([join(sh, "{}.ccp4".format(probe)) for sh in shrunk_hots])
            ge = GridEnsemble(ge_paths)
            ge.get_ensemble_array()
            save_path = join(self.root_dir, "{}_{}.p".format(self.ensemble_ID, probe))
            ge.save_gridensemble(save_path)

    def save_ensemble_data(self):
        pass


class MultipleEnsemble(Ensemble):
    """
    Class that aligns and performs comparisons on multiple EnsembleResults (e.g. target protein and >=1 off-targets in same family.
    """
    def __init__(self, root_dir, ref_id, ensemble_dict, align_references=True, references_rmsds=None):
        """
        
        :param str root_dir: path to the top-level directory where the SIENA output and EnsembleResults ouputs will be stored
        :param pandas.Series ref_id: a row from a the protein_data of the reference ResultEnsemble
        :param ensemble_dict: dictionary in the form {EnsembleResult.ensemble_ID : <EnsembleResult instance>}
        :param bool align_references: True if ensembles have not been previously aligned. False if supplying pre-aligned ensembles
        :param list references_rmsds: Optional. If the ensembles have been pre-aligned, can be supplied for extra info. If align_refereces
                                       is False, this will get updated with the values from the alignment performed by this class.
        """
        Ensemble.__init__(self, root_dir, ref_id)
        self.ensembles = ensemble_dict
        self.reference_ensemble = ensemble_dict[self.reference_ID["Ensemble ID"]]
        self.aligned_references = align_references
        self.ensemble_references_RMSD = references_rmsds
        self.savedir = join(self.root_dir, "EnsembleComparisons")

    def align_references(self):
        """
        For all the EnsembleResults in the MultipleEnsemble, aligns the reference structure of each EnsembleResult to that
        of the reference result of the MultipleEnsemble
        :return: 
        """
        # Get the binding sites of the reference structure for each ensemble
        self.reference_ensemble.reference_binding_site = self.get_binding_site(self.reference_ID)

        for key, ens in self.ensembles.items():
            print(key, ens.reference_ID)
            ens.reference_binding_site = self.get_binding_site(ens.reference_ID)

        ref_bs = self.reference_ensemble.reference_binding_site
        self.reference_ensemble._save_protein(ref_bs)

        rmsds = []
        for key, ens in self.ensembles.items():
            chain_superposition = Protein.ChainSuperposition(self.alignment_settings)
            (rmsd, transformation) = chain_superposition.superpose(ref_bs.protein.chains[0],
                                                                   ens.reference_binding_site.protein.chains[0],
                                                                   ref_bs)

            #ens.process_ensemble_proteins()

            rmsds.append("RMSD ref: {} {} query: {} {} = {:6.4f}".format(self.reference_ensemble.ensemble_ID,
                                                                         self.reference_ensemble.reference_ID["ID"],
                                                                         ens.ensemble_ID,
                                                                         ens.reference_ID["ID"],
                                                                         rmsd))
        print(rmsds)
        self.aligned_references = True
        self.ensemble_reference_RMSD = rmsds

    def binding_site_ensembles(self):
        """
        For each ensemble in EnsembleResult
        :return: 
        """
        # Use the fact that the EnsembleResult.run_hotspots() method returns the paths to all the hotspots in the ensemble.
        # It will not recalculate hotspots if they are already present (Luigi will recognise them as localtargets).
        hot_path_dic = {ens.ensemble_ID: ens.run_hotspots() for ens in self.ensembles.values()}

        # The shrink_hotspots() function is called from this class (not the individual EnsembleResults), so the shrunk hotspots of
        # all the ensembles will be same size (set by the boundaries of the reference structure of the MultipleEnsemble).
        shrunk_hot_paths = {key: self.shrink_hotspots(hot_paths) for key, hot_paths in hot_path_dic.items()}

        return shrunk_hot_paths



















from __future__ import print_function, division
import pandas as pd
from os import mkdir
from os.path import join, exists, basename, dirname
from glob import glob
from ccdc.io import MoleculeReader, MoleculeWriter
from ccdc.protein import Protein
from Ensemble import Ensemble, EnsembleResult
import luigi
from run_hotspots import ParalleliselRunner
import collections
import numpy as np



class SIENAreader:
    """
    A class to read SIENA output (SIENA taken from https://proteins.plus/) and perform basic filtering.
    """
    def __init__(self, siena_dir, ens_name, ref_pdb, out_dir=""):
        """
        Given the directory storing the SIENA output, read it in and perform basic filtering
        :param str siena_dir:  path to where SIENA output is stored.
        """
        self.out_dir = out_dir
        self.direc = siena_dir
        self.ensemble_name = ens_name
        self.reference_pdb = ref_pdb
        self.reference_structure = None
        # self.lig_dir = join(self.direc, "ligands")
        # self.pdb_dir = join(self.direc, "pdb_files")
        self.lig_dir = join(self.direc, "ligand")
        self.pdb_dir = join(self.direc, "ensemble")
        #self.result_stats = pd.read_csv(join(self.direc,"result_table", "resultStatistic.csv"), sep=";")
        self.result_stats = pd.read_csv(join(self.direc, "resultStatistic_fragments.csv"), sep=";")
        self.alignment_file = None
        self.proteins_df = None
        self.aa_dic = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS", "I": "ILE", "K" : "LYS",
                       "L": "LEU", "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
                       "W": "TRP", "Y": "TYR"}


    def read_alignment(self):
        """
        Reformats the SIENA output file in a way that can be easily read by pandas. Outputs the alignment as a dataframe.
        :return: pandas.dataframe object
        """
        # Want to avoid modifying the original output file, so we create a copy
        new_path = join(self.direc, "pd_readable_alignment.csv")
        self.alignment_file = new_path

        with open(join(self.direc, "alignment", "alignment.txt"), "r") as input:
            with open(self.alignment_file, "w") as output:
                l = input.readlines()
                # The first line says "Alignment:" - we want to get rid of it
                l_out = l[1:]
                [output.write(l_o.replace(" ", "")) for l_o in l_out]

        ali_df = pd.read_csv(self.alignment_file,  sep="|")
        # removes the "Nan" column created at the end of the dataframe by the final separator.
        ali_df = ali_df.iloc[:, :-1]

        return ali_df

    def _get_aa_from_alignment(self, s):
        """
        Each cell of the alignment_df is formatted "<1 letter amino acid code>-<residue number>-<chain>"
        This function returns the amino acid, converted to a 3-letter code (for compatibility with CSD API).
        :param str s: contents of cell in the alignment dataframe. If not correct type of cell, returns nothing.
        :return: str or None
        """
        # Checks that the string is in the format described above.
        if len(s.split("-")) > 2:
            # get the 1 letter aa code from the string
            aa_letter = s.split("-")[0].strip()
            # Convert to 3-letter aa code
            aa_full = self.aa_dic[aa_letter]
            return aa_full
        else:
            return

    def get_chains(self, row):
        """
        The SIENA output also has a "chains" section, but this way we get the chains involved in the binding site only
        :param row: a row in the alignment dataframe
        :return: list
        """
        # row.values[0] is the PDB code
        #vals = row.values[1:]
        vals = row[1:]
        print(vals)
        # get a list of chains for each amino acid in the binding site
        allchains = [s.split("-")[2].strip() for s in vals]
        chains = list(set(allchains))
        return chains

    def get_binding_site_residues(self, row):
        """
        Converts the binding site information in the SIENA output to a format that can be passed to the CSD API
        :param row: a row in the alignment dataframe
        :return: 
        """
        # row.values[0] is the PDB code
        vals = row[1:]
        # get the amino acid position in the structure
        pos = [s.split("-")[1].strip() for s in vals]
        # get the amino acids and convert to 3-letter code
        aas = [self._get_aa_from_alignment(s) for s in vals]
        # get the chain for each amino acid
        chains = [s.split("-")[2].strip() for s in vals]
        # format
        bs_res = ["{0}:{1}{2}".format(chains[i], aas[i], pos[i]) for i in range(len(vals))]

        return bs_res

    def make_row_dic(self, idx, row):
        """
        Compiles each row of the dataframe that gets passed to the EnsembleResult.
        :param int idx: row index
        :param pandas.Series row: 
        :return: dictionary
        """
        try:
            chains = self.get_chains(row)
            row_dic = {"Filename": join(self.pdb_dir, "{}_{}.pdb".format(row[0], idx+1)),
                       # Assume that the combination of PDB ID and binding site residues is unqiue. Possibly problematic for NUDT5?
                       "ID": "{}-{}".format(row[0], "".join(chains)),
                       "PDB ID": row[0],
                       "Chains": chains,
                       "Binding site": self.get_binding_site_residues(row)}
            return row_dic
        except AttributeError:
            return

    def get_ensemble_protein_df(self):
        """
        Creates the dataframe that holds the ensemble information and will be propagated to the EnsembleResult.
        :return: 
        """
        ali_df = self.read_alignment()
        #ali_df = ali_df.iloc[:, :-1]

        # Create a dataframe with the column names needed by EnsembleResult
        new_rows = ["ID", "Filename", "PDB ID", "Chains", "Binding site"]
        a = pd.DataFrame(columns=new_rows, dtype="object")

        ali_df_vals = ali_df.values
        # Fill dataframe
        for idx, r in enumerate(ali_df_vals):
            print(idx, r)
            a.loc[len(a)] = self.make_row_dic(idx, r)
        # All structures should have the same ensemble ID (e.g. target name)
        a["Ensemble ID"] = self.ensemble_name

        # Get the information on the reference structure to pass on to EnsembleResult
        ref_row = a[a["PDB ID"]==self.reference_pdb].iloc[0]
        #ref_pdb = self.reference_pdb
        print(ref_row["Chains"])
        ref_chains = "".join(ref_row["Chains"])
        print(ref_row, self.reference_pdb, ref_chains)
        a.to_csv(join(self.out_dir, "ensemble_data.csv"))
        ref_id = a.loc[a["ID"] == "{}-{}".format(self.reference_pdb, "".join(ref_chains))]
        self.reference_structure = ref_id

        return a

    def find_largest_ligand(self):
        """
        Looks for the largest ligand returned in the SIENA ensemble.
        :return: 
        """
        # Get the ligands for the proteins returned by SIENA
        print(self.lig_dir)
        mol_paths = glob(join(self.lig_dir, "*.sdf"))
        print(mol_paths)
        mols = MoleculeReader(mol_paths)

        # Get a dictionary of the molecule_ID and the filename
        mw_dict = {basename(m_fname): m.molecular_weight for m_fname, m in zip(mols.file_name, mols)}

        print(mw_dict)
        # Get the filename of the largest ligand:
        try:
            largest_lig = sorted(((value,key) for (key,value) in mw_dict.items()), reverse=True)[0][1]
            print(largest_lig)
            return largest_lig

        except IndexError:
            print("SIENA found no ligands for ensemble {}".format(self.ensemble_name))
            return


    def get_largest_binding_site(self):
        """
        Returns the binding site created within 6.5A of the largest ligand
        :return: 
        """
        lig_fname = self.find_largest_ligand()
        lig = MoleculeReader(join(self.lig_dir, lig_fname))[0]
        prot = Protein.from_file(join(self.pdb_dir, lig_fname.replace("sdf", "pdb")))
        bs = Protein.BindingSiteFromMolecule(protein=prot, molecule=lig, distance=6.5)
        return bs

    def output_ensemble(self):
        print("outputting ensemble", self.ensemble_name)
        prot_df = self.get_ensemble_protein_df()
        e = EnsembleResult(root_dir=self.out_dir,
                           ref_id=self.reference_structure.squeeze().to_dict(),
                           df=prot_df)
        e.ensemble_ID = self.ensemble_name
        if not exists(self.out_dir):
            mkdir(self.out_dir)
        return e


class KLIFSReader:
    """
    Class that reads in the KLIFS downloaded data.
    """
    def __init__(self, ensemble_directory, ens_name):
        self.root_dir = ensemble_directory
        self.ensemble_name = ens_name
        self.ensemble_data = pd.read_csv(join(self.root_dir, "overview.csv"))
        # if exists(join(self.root_dir, "subset_overview.csv")):
        #     self.ensemble_data = pd.read_csv(join(self.root_dir, "subset_overview.csv"))
        # else:
        #     self.ensemble_data = pd.read_csv(join(self.root_dir, "overview.csv"))

    def find_largest_ligand(self):
        """
        Finds the ligand with the highest molecular weight in the ensmeble.
        :return: 
        """
        lig_paths = glob(join(self.root_dir, '*', 'ligand.mol2'))
        ligs = MoleculeReader(lig_paths)
        mw_dic = {m_fname: m.molecular_weight for m_fname, m in zip(ligs.file_name, ligs)}
        largest_lig = sorted(((value, key) for (key, value) in mw_dic.items()), reverse=True)[0][1]
        print('Largest ligand for {}: {}'.format(self.ensemble_name, mw_dic[largest_lig]))
        # Note - this is only the path to the largest_lig
        return largest_lig

    def get_protein_paths(self, mode='protein'):
        prot_paths = []
        for i, row in self.ensemble_data.iterrows():
            print(i, row)
            if row['alt'] != ' ':
                stri = "{}_alt{}_chain{}".format(row['pdb'], row['alt'], row['chain'])
                #only take alt A
                #stri = "{}_altA_chain{}".format(row['pdb'], row['chain'])
            else:
                stri = "{}_chain{}".format(row['pdb'], row['chain'])
            prot_paths.append(join(self.root_dir, stri, '{}.mol2'.format(mode)))

        return prot_paths


    def save_pdbs(self, paths_list = None, data_frame=None):
        """
        Takes the KLIFS output and saves as pdbs
        :return: 
        """
        if not paths_list:
            complexes = self.get_protein_paths( mode='complex')
        else:
            complexes = paths_list
        comp_dir = join(self.root_dir, 'Complexes')
        if not exists(comp_dir):
            mkdir(comp_dir)
        print('reading in proteins....')
        for f_name in complexes:
            prot = Protein.from_file(f_name)
            prot.detect_ligand_bonds()
            new_name = join(comp_dir, "{}.pdb".format(basename(dirname(f_name))))
            with MoleculeWriter(new_name) as protein_writer:
                protein_writer.write(prot)

    def get_subset(self):
        """
        Returns a subset of paths based on what is 
        :return: 
        """
        df = self.ensemble_data
        # First, get the most common sequence (we assume that's the true one).
        seq_counter = collections.Counter(df['pocket'])
        most_common_list = seq_counter.most_common()
        most_common_seq = most_common_list[0][0]

        # If there is more than one sequence, say what it is and truncate the dataframe
        if len(most_common_list) >1:
            print("Chosen binding site sequence: \n {} \n Next most commonly occurring: \n {}".format(most_common_list[0],
                                                                                                        most_common_list[1]))
            df = df[df['pocket'] == most_common_seq]

        # Remove structures in complex with the same ligands. Keep only the highest resolution one.
        ligs = list(set(df['orthosteric_PDB']))
        idx_list = []
        for lig in ligs:
            sli = df[df['orthosteric_PDB'] == lig]
            # get the entry with the highest resolution, or if there is a tie, take the first.
            idx = sli[sli['resolution']==sli['resolution'].min()].index[0]
            idx_list.append(idx)
        print(idx_list[:60], len(idx_list[:60]))
        self.ensemble_data = self.ensemble_data.iloc[idx_list[:60]] # Get the first 60 structures from all
        #self.ensemble_data.to_csv(join(self.root_dir, 'subset_overview.csv'))

    def no_allosteric(self):
        """
        Removes all allosteric ligands
        :return: 
        """
        df = self.ensemble_data

        # Make a separate directory for the output:
        new_dir = join(self.root_dir, "no_allosteric")
        if not exists(new_dir):
            mkdir(new_dir)
        prot_paths = self.get_protein_paths(mode="complex")
        df['ID'] = [dirname(p) for p in prot_paths]


        hot_paths = {}
        for path in df["ID"]:
            tar_path = join(path, "fullsize_hotspots_100000", "binding_site_maps", "out")
            if exists(tar_path):
                hot_paths[path] = tar_path

        no_allo_hot_paths = {}
        for h in hot_paths.keys():
            row =  df.loc[df['ID']==h].squeeze()
            print(row)
            if row['allosteric_PDB'] == "-":
                no_allo_hot_paths[h] = hot_paths[h]

        from GridEnsemble import GridEnsemble

        probes = ["donor", "acceptor", "apolar"]
        for probe in probes:
            probe_paths = [join(path, "{}.ccp4".format(probe)) for path in no_allo_hot_paths.values()]
            print(probe_paths)
            ge = GridEnsemble(probe_paths)
            ge.get_ensemble_array()
            save_path = join(new_dir, "{}_{}.p".format(self.ensemble_name, probe))
            ge.save_gridensemble(save_path)

        s_paths = [join(p, "complex.mol2") for p in no_allo_hot_paths.keys()]
        self.save_pdbs(s_paths)

        return

    def get_ensemble(self, nrotations, charged=False):
        largest_lig = self.find_largest_ligand()
        lig = MoleculeReader(largest_lig)[0]
        prot = Protein.from_file(join(dirname(largest_lig), 'protein.mol2'))
        bs = Protein.BindingSiteFromMolecule(protein=prot,
                                             molecule=lig,
                                             distance=6.5)

        # prot_paths = glob(join(self.root_dir, '*', 'protein.mol2'))
        prot_paths = self.get_protein_paths()
        print(prot_paths)
        print(self.ensemble_name, len(prot_paths))
        luigi.build([ParalleliselRunner(prot_paths, nrotations, charged, data_source='KLIFS')], local_scheduler=True, workers=30)
        #luigi.build([ParalleliselRunner(prot_paths, nrotations, charged)], local_scheduler=True,
                    #workers=30)
        hot_paths = [join(dirname(in_pdb), "fullsize_hotspots_{}".format(nrotations), "out.zip") for in_pdb in prot_paths]
        return hot_paths


if __name__ == "__main__":
    import tempfile
    import json
    from GridEnsemble import GridEnsemble

    tempfile.tempdir = "/home/jin76872/Desktop/Mih/Data/tmp_superstar_ghecom"
    #ens = ['p38a', 'PIM1', 'Erk2']
    ens = ['Erk2']
    main_dir = "/home/jin76872/Desktop/Mih/Data/ACS_fall_meeting_2019_slides"
    nrot = 100000

    for e in ens:
        kr = KLIFSReader(ensemble_directory=join(main_dir, e),
                         ens_name=e)
        #kr.no_allosteric()
        kr.get_subset()
        #kr.save_pdbs()

        #kr.find_largest_ligand()

        hot_paths = kr.get_ensemble(nrotations=nrot, charged=False)


        ref_kr = KLIFSReader(ensemble_directory=join(main_dir, 'CK2a1'), ens_name= 'CK2a1')
    # # CK2_kr.save_pdbs()
    #
    #
        largest_ligand = ref_kr.find_largest_ligand()
        prot = Protein.from_file(join(dirname(largest_ligand), 'protein.mol2'))
        bs = Protein.BindingSiteFromMolecule(protein=prot,
                                             molecule=MoleculeReader(largest_ligand)[0],
                                             distance=7.0)

        s_paths_file = join(main_dir, "shrunk_hot_paths_{}.json".format(nrot))
        if exists(s_paths_file):
            with open(s_paths_file, "r") as f:
                s_paths_dict = json.load(f)
        else:
            s_paths_dict = {}

        ensemble = Ensemble(root_dir=join(main_dir, e))
        ensemble.reference_binding_site = bs
        #hot_paths = glob(join(ensemble.root_dir, '*', "fullsize_hotspots_{}".format(nrot), "out.zip"))
        s_paths = ensemble.shrink_hotspots(hotspot_paths=hot_paths,
                                 padding=2.0)
        s_paths_dict[e] = s_paths
        for probe in ["donor", "acceptor", "apolar"]:
            paths = [join(t, '{}.ccp4'.format(probe)) for t in s_paths]
            gr = GridEnsemble(paths)
            gr.get_ensemble_array()
            gr.save_gridensemble(join(main_dir, e, '{}_{}.p'.format(probe, nrot)))

    spd = json.dumps(s_paths_dict, sort_keys=True, indent=4, separators=(',', ': '))
    with open(join(main_dir, "shrunk_hot_paths_{}.json".format(nrot)), "w") as f:
        f.write(spd)

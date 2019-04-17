from __future__ import print_function, division
import pandas as pd
from os.path import join
from Ensemble import EnsembleResult


class SIENAreader:
    """
    A class to read SIENA output (SIENA taken from https://proteins.plus/) and perform basic filtering.
    """
    def __init__(self, siena_dir, ens_name, out_dir=""):
        """
        Given the directory storing the SIENA output, read it in and perform basic filtering
        :param str siena_dir:  path to where SIENA output is stored.
        """
        self.out_dir = out_dir
        self.direc = siena_dir
        self.ensemble_name = ens_name
        self.reference_structure = None
        self.lig_dir = join(self.direc, "ligand")
        self.pdb_dir = join(self.direc, "ensemble")
        self.result_stats = pd.read_csv(join(self.direc, "resultStatistic.csv"), sep=";")
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

        with open(join(self.direc, "alignment.txt"), "r") as input:
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
        vals = row.values[1:]
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
        vals = row.values[1:]
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
        chains = self.get_chains(row)
        row_dic = {"Filename": join(self.pdb_dir, "{}-{}.pdb".format(row[0], idx+1)),
                   # Assume that the combination of PDB ID and binding site residues is unqiue. Possibly problematic for NUDT5?
                   "ID": "{}-{}".format(row[0], "".join(chains)),
                   "PDB ID": row[0],
                   "Chains": chains,
                   "Binding site": self.get_binding_site_residues(row)}
        return row_dic

    def get_ensemble_protein_df(self):
        """
        Creates the dataframe that holds the ensemble information and will be propagated to the EnsembleResult.
        :return: 
        """
        ali_df = self.read_alignment()
        ali_df = ali_df.iloc[:, :-1]

        # Create a dataframe with the column names needed by EnsembleResult
        new_rows = ["ID", "Filename", "PDB ID", "Chains", "Binding site"]
        a = pd.DataFrame(columns=new_rows, dtype="object")

        # Fill dataframe
        for idx, r in ali_df.iterrows():
            a.loc[len(a)] = self.make_row_dic(idx, r)
        # All structures should have the same ensemble ID (e.g. target name)
        a["Ensemble ID"] = self.ensemble_name

        # Get the information on the reference structure to pass on to EnsembleResult
        ref_row = ali_df.columns
        ref_pdb = ref_row.values[0]
        ref_chains = self.get_chains(ref_row)

        ref_id = a.loc[a["ID"] == "{}-{}".format(ref_pdb, "".join(ref_chains))]
        self.reference_structure = ref_id

        return a

    def output_ensemble(self):
        e = EnsembleResult(root_dir=self.out_dir,
                           ref_id=self.reference_structure.squeeze(),
                           df=self.get_ensemble_protein_df())
        e.ensemble_ID = self.ensemble_name
        return e
import pandas as pd
import collections
from pathlib import Path


def get_protein_paths(df, stem_directory, mode='protein'):
    prot_paths = []
    for i, row in df.iterrows():
        print(i, row)
        if row['alt'] != ' ':
            stri = "{}_alt{}_chain{}".format(row['pdb'], row['alt'], row['chain'])
            # only take alt A
            # stri = "{}_altA_chain{}".format(row['pdb'], row['chain'])
        else:
            stri = "{}_chain{}".format(row['pdb'], row['chain'])
        prot_paths.append(Path(stem_directory, stri, '{}.mol2'.format(mode)))

    return prot_paths

def get_subset(summary_csv, new_df_path):
    """
    Returns a subset of paths based on what is 
    :return: 
    """
    df = pd.read_csv(summary_csv)


    # Remove structures in complex with the same ligands. Keep only the highest resolution one.
    ligs = list(set(df['orthosteric_PDB']))
    idx_list = []
    for lig in ligs:
        sli = df[df['orthosteric_PDB'] == lig]
        # get the entry with the highest resolution, or if there is a tie, take the first.
        idx = sli[sli['resolution'] == sli['resolution'].min()].index[0]
        # Don't include allosteric ligands

        idx_list.append(idx)
    print(idx_list)
    df = df.iloc[idx_list]

    # First, get the most common sequence (we assume that's the true one).
    seq_counter = collections.Counter(df['pocket'].values)
    most_common_list = seq_counter.most_common()
    most_common_seq = most_common_list[0][0]

    # If there is more than one sequence, say what it is and truncate the dataframe
    if len(most_common_list) > 1:
        print("Chosen binding site sequence: \n {} \n Next most commonly occurring: \n {}".format(most_common_list[0],
                                                                                                  most_common_list[1]))
        df = df[df['pocket'] == most_common_seq]

    # Get rid of any leftover allosteric ligands?
    prot_paths = get_protein_paths()

    df.to_csv(new_df_path, index=False)

    return df

if __name__ == "__main__":
    cdk2_csv = Path("/home/jin76872/Desktop/Mih/Data/ensemble_maps_validation/CDK2/overview.csv")
    out_path = Path("/home/jin76872/Desktop/Mih/Data/ensemble_maps_validation/CDK2/removed_duplicates.csv")
    df_new = get_subset(cdk2_csv, out_path)
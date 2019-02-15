from ccdc.protein import Protein
import os
from os.path import join, basename, exists
from glob import glob
from hotspots.calculation import Runner
from hotspots import hs_io
#from ccdc.cavity import Cavity
import sys


def make_savedir(stem, identifier):
    save_dir = join(stem, identifier, "{}_hotspot_result".format(identifier))
    if not exists(save_dir):
        os.mkdir(save_dir)
    return save_dir


def prepare_protein(prot_path, ch=["A", "B"]):
    prot = Protein.from_file(prot_path)
    for chain in prot.chains:
        print(chain.identifier)
        if chain.identifier not in ch:
            prot.remove_chain(chain.identifier)
    prot.remove_all_waters()
    for ligand in prot.ligands:
        prot.remove_ligand(ligand.identifier)
    prot.remove_all_metals()
    prot.add_hydrogens()
    return prot


def calc_hotspot(path, prot_name, method, nrot=3000):
    """
    param: path str, path to prepared protein
    """
    protein = prepare_protein(path)
    h = Runner()

    settings = h.Settings()
    settings.nrotations = nrot
    settings.apolar_translation_threshold = 15
    settings.polar_translation_threshold = 15
    settings.sphere_maps=False

    result = h.from_protein(protein=protein,
                            charged_probes=False,
                            probe_size=7,
                            buriedness_method=method,
                            cavities=None,
                            nprocesses=3,
                            settings=settings)
    #out = make_savedir(prot_name)
    out = os.getcwd()
    with hs_io.HotspotWriter(out, visualisation="pymol", grid_extension=".ccp4", zip_results=False) as writer:
        writer.write(result)
    return result

#p= "/media/mihaela/BIG_DATA_IOT_BLOCKCHAIN/PhD/SABS_year2/ensemble_example/BRD1A-x050_event1_aligned_mal_mega.pdb"

p = sys.argv[1]
prot_name = basename(p).split(".")[0]
res = calc_hotspot(p, prot_name, "ghecom")
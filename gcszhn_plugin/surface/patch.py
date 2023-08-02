import numpy as np

from pymol import cmd
from scipy.spatial.distance import cdist
from .mesh_utils import Mesh

__all__ = ['extract_patch', 'patch_seq']

# 三字母氨基酸缩写和单字母缩写的对应关系，字典
# https://www.bioinformatics.org/sms/iupac.html
amino_code = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y"
}



def extract_patch(patch_name:str, patch_ply_name:str, model_name: str = "(all)", remove_model: bool = False,
                  distance_threshold: float = 4.0):
    atoms_coords = np.array([atom.coord for atom in cmd.get_model(model_name).atom])
    atoms_ids = np.array([atom.id for atom in cmd.get_model(model_name).atom])
    mesh = Mesh.create_mesh()
    mesh.load_mesh(patch_ply_name)
    atoms_dists = cdist(atoms_coords, mesh.vertices).min(axis=1)
    atoms_selected = atoms_ids[atoms_dists < distance_threshold]
    ids_selected = ",".join([str(i) for i in atoms_selected])
    if model_name == "(all)":
        cmd.select("selected_atoms", f"id " + ids_selected)
    else:
        cmd.select("selected_atoms", f"model {model_name} and id {ids_selected}")
    cmd.extract(patch_name, "selected_atoms")
    if remove_model:
        cmd.delete(model_name)
    cmd.delete("selected_atoms")


def patch_seq(patch_name):
    resi_set = set()
    seq_list = []
    for atom in cmd.get_model(patch_name).atom:
        resn = atom.resn
        resi = int(atom.resi)
        if resi in resi_set:
            continue
        resi_set.add(resi)
        while resi > len(seq_list) + 1:
            seq_list.append("-")
        seq_list.append(amino_code[resn])
    return "".join(seq_list)


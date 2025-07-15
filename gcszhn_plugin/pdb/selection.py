from pymol import cmd
from ..utils import register_pymol_cmd


@register_pymol_cmd
def findseq(seq: str, name: str = 'sele', obj: str = 'all'):
    """
    Find a sub sequence in the PDB object.
    The sequence is a string of one-letter amino acid code.
    """
    aa_one2three = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
        'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
        'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
        'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }
    seq = seq.upper()
    chain_ids = cmd.get_chains(obj)
    for chain_id in chain_ids:
        chain = cmd.get_model(f'{obj} and chain {chain_id}')
        res_first_atoms = [chain.atom[r[0]] for r in chain.get_residues()]
        s, e = 0, 0
        while e < len(res_first_atoms):
            if res_first_atoms[e].resn == aa_one2three[seq[e-s]]:
                e += 1
                if e - s == len(seq):
                    cmd.select(name, f'{obj} and chain {chain_id} and resi {res_first_atoms[s].resi}-{res_first_atoms[e-1].resi}')
                    return
            else:
                s += 1
                e = s
    else:
        print(f'Cannot find sequence {seq} in {obj}')
        return
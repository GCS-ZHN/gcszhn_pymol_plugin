import pandas as pd

from pymol import cmd
from ..utils import register_pymol_cmd, int_array_to_str

__all__ = ['copy_selection']


@register_pymol_cmd
def copy_selection(selection: str = "(all)", mode="tab"):
    """
    Copy selection to clipboard.

    Parameters
    ----------
    selection : str, optional
        Selection.
    mode : str, optional
        Copy mode. One of "tab", "rf_hotspot", "range".
        tab: copy chain, resi, resn to clipboard in TSV format.
        rf_hotspot: copy chain and resi in RFdiffusion hotspot format.
        range: copy chain and resi in numerial range format.
    """
    def _gen():
        for atom in cmd.get_model(selection).atom:
            yield atom.chain, atom.resi, atom.resn
    
    data = pd.DataFrame(_gen(), columns=['chain', 'resi', 'resn'])
    # remove duplicated atoms
    data = data.drop_duplicates(subset=['chain', 'resi', 'resn'])
    if mode == "tab":
        data.to_clipboard(index=False)
    elif mode == "rf_hotspot":
        resi = data['chain'] + data['resi']
        resi = pd.DataFrame(resi, columns=['resi'])
        resi.T.to_clipboard(index=False, sep=",", header=False)
    elif mode == "range":
        data_group = data.groupby("chain")
        resi_range = data_group.apply(lambda d: int_array_to_str(d['resi'].astype(int), d.name))
        resi_range.to_clipboard(index=False, header=False)
    else:
        raise ValueError(f"Unknown mode: {mode}")


@register_pymol_cmd
def split_by_chain(obj_name: str = None):
    """
    Split object into multiple objects by chain.

    Parameters
    ----------
    obj_name : str, optional
        Object name. If None, all objects will be processed.
    """
    object_names = cmd.get_object_list()
    if obj_name and obj_name in object_names:
        object_names = [obj_name]
    for obj in object_names:
        chains = cmd.get_chains(obj)
        if len(chains) > 1:
            for chain in chains:
                new_obj_name = f"{obj}_{chain}"
                cmd.create(new_obj_name, f"{obj} and chain {chain}")
            cmd.delete(obj)

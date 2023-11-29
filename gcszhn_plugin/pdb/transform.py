import pandas as pd

from pymol import cmd
from ..utils import register_pymol_cmd, int_array_to_str

__all__ = ['copy_selection']


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

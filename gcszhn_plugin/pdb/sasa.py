from typing import Dict
from pymol import cmd
from ..utils import local_setting, residue_format, register_pymol_cmd, residue_with_CA

__all__ = ["set_sasa_color", "get_sasa"]


@register_pymol_cmd
def get_sasa(
        selection: str,
        solvent_radius: float = 1.4,
        load_b: int = 0,
        dot_density: int = 2) -> float:
    """
    Calculate solvent accessible surface area (SASA)
    for specific selection.
    """
    with local_setting(
        dot_solvent='on', 
        dot_density=dot_density, 
        solvent_radius=solvent_radius):
        return cmd.get_area(selection, load_b=load_b)


@register_pymol_cmd
def set_sasa_color(
        selection:str = "(all)", 
        level: str = "A", 
        palette: str = 'red_white_blue',
        minimum: float = None,
        maximum: float = None):
    """
    Annotated color by sasa.

    Parameters
    ---------------
    - selection: str
        pymol selection str.
    - level: str
        'A' for atomic level,
        'R' for residual level,
        'C' for chain level.
    - palette: str
        color palette of pymol.
    
    """
    
    if type(minimum) is not type(maximum):
        raise ValueError("Please specific minimum and maximum both or not!")
    
    auto_bound = minimum is None
    
    if auto_bound:
        minimum = float('inf')
        maximum = -minimum

    sasa_buffer = dict()

    get_sasa(selection, load_b=1)

    def set_sasa(resi, chain, b):
        if level == "R":
            res_sel_str = f'resi {resi} and chain {chain} and {selection}'
        elif level == "C":
            res_sel_str = f'chain {chain} and {selection}'

        if level == 'A':
            v = b
        else:
            if res_sel_str not in sasa_buffer:
                sasa_buffer[res_sel_str] = sum(atom.b for atom in cmd.get_model(res_sel_str).atom)
            v = sasa_buffer[res_sel_str]
        
        if auto_bound:
            nonlocal minimum, maximum
            minimum = min(minimum, v)
            maximum = max(maximum, v)
        return v
    
    cmd.alter(
        selection,
        f'b = set_sasa(resi, chain, b)',
        space={'set_sasa': set_sasa})
    cmd.spectrum(
        'b',
        palette=palette,
        selection=selection,
        minimum=minimum,
        maximum=maximum)


def get_sasa_by_res(selection: str) -> Dict[str, float]:
    get_sasa(selection=selection, load_b=1)
    selection = residue_with_CA(selection)
    sasa_dict = dict()

    def _update_res_sasa(atom):
        residue = residue_format(atom.resn, atom.resi, atom.chain, selection)
        sasa_dict[residue] = sum(atom.b for atom in cmd.get_model(residue).atom)
    
    cmd.iterate(selection, _update_res_sasa)
    return sasa_dict

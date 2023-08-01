import numpy as np
import pandas as pd

from pymol import cmd
from colour import Color
from importlib import resources

# load scale.csv under current module dir
with resources.open_text(__package__, "hydro_scale.csv") as f:
    HYDRO_SCALE_MAP = pd.read_csv(f, index_col=0)


def get_hydro(resn, scale_name='Ja', ignore_miss_res: bool = True):
    scale_map = HYDRO_SCALE_MAP[scale_name]
    if resn in scale_map.index:
        return scale_map[resn]
    else:
        if not ignore_miss_res:
            raise KeyError(f"Undefined resn {resn} in hydrophobicity scale.")
        return 0


def avail_hydro_scales():
    """
    List available hydrophobicity scales.
    """
    for scale_name in HYDRO_SCALE_MAP.columns:
        print(scale_name)


def colormap(
        val: float,
        minimum: float = 0.0,
        maximum: float = 1.0,
        min_color: str = 'blue',
        max_color: str = 'red',
        step: int = 100) -> Color:
    if not isinstance(step, int) or step <= 0:
        raise TypeError('step must be positive integer')
    min_color = Color(min_color)
    max_color = Color(max_color)
    val_range = np.arange(minimum, maximum, (maximum - minimum) / step)
    color_range = min_color.range_to(max_color, 100)
    for val_up_bound, val_color in zip(val_range, color_range):
        if val <= val_up_bound:
            return val_color
    return val_color


def set_hydro_color_v2(
        selection: str = '(all)',
        palette: str = 'blue_white_red',
        scale_name: str = 'Ja',
        minimum: float = None,
        maximum: float = None,
        level: str = 'A',
        with_sasa: bool = True):
    """
    Annotate hydrophobicity color of each
    residue. this function implmented by
    `pymol.cmd.alter` and `pymol.cmd.spectrum`.

    Parameters
    ----------
    selection: str
        pymol selection string.
    scale_name: str
        Name of hydrophobicity scale.
        list all by `available_hydro_scales()`.
    palette: str
        interal pymol color palette.
        more available palette is available
        by `help spectrum`.
    """
    if type(minimum) is not type(maximum):
        raise ValueError("Please specific minimum and maximum both or not!")
    
    auto_bound = minimum is None
    
    if auto_bound:
        minimum = float('inf')
        maximum = -minimum

    dot_solvent = cmd.get('dot_solvent')
    sasa_buffer = dict()

    if with_sasa:
        cmd.set('dot_solvent', 'on')
        cmd.get_area(selection, load_b=1)

    def set_hydro(resn, resi, chain, b, scale_name='Ja'):
        nonlocal minimum, maximum
        v = get_hydro(resn, scale_name=scale_name)
        if with_sasa:
            if level == 'R':
                res_sel_str = f'resi {resi} and chain {chain} and {selection}'
            elif level == 'C':
                res_sel_str = f'chain {chain} and {selection}'
            
            if level == 'A':
                sasa = b
            else:
                if res_sel_str not in sasa_buffer:
                    sasa_buffer[res_sel_str] = sum(atom.b for atom in cmd.get_model(res_sel_str).atom)
                sasa = sasa_buffer[res_sel_str]
            v *= sasa
        
        if auto_bound:
            minimum = min(minimum, v)
            maximum = max(maximum, v)
        return v
    
    cmd.alter(
        selection,
        f'b = set_hydro(resn, resi, chain, b, scale_name="{scale_name}")',
        space={'set_hydro': set_hydro})
    cmd.spectrum(
        'b',
        palette=palette,
        selection=selection,
        minimum=minimum,
        maximum=maximum)

    if with_sasa:
        cmd.set('dot_solvent', dot_solvent)


def set_hydro_color_v1(
        selection: str = '(all)',
        scale_name: str = 'Ja',
        min_color: str = 'blue',
        max_color: str = 'red'):
    """
    Annotate hydrophobicity color of each
    residue.

    Parameters
    -----------
    selection: str
        pymol selection string.
    scale_name: str
        Name of hydrophobicity scale.
        list all by `available_hydro_scales()`.
    min_color: str
        Color of minimum hydrophobicity.
    max_color: str
        Color of maximum hydrophobicity.
    """
    scale_map = HYDRO_SCALE_MAP[scale_name]
    minimum = scale_map.min()
    maximum = scale_map.max()
    for resn, score in scale_map.items():
        res_color = colormap(
            score,
            minimum=minimum,
            maximum=maximum,
            min_color=min_color,
            max_color=max_color)
        cmd.select('resn_sel', f'resn {resn} and {selection}')
        colorn = f'{resn}_color_sel_{selection}'
        cmd.set_color(
            colorn,
            res_color.get_rgb())
        cmd.color(colorn, 'resn_sel')
        cmd.delete('resn_sel')

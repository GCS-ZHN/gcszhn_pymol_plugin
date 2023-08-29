import numpy as np
import pandas as pd

from pymol import cmd
from colour import Color
from importlib import resources
from .sasa import get_sasa_by_res
from ..utils import register_pymol_cmd, residue_with_CA, residue_format


__all__ = [
    "set_hydro_color_v2",
    "set_hydro_color_v1",
    "avail_hydro_scales",
    "set_hydration_color"]

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


@register_pymol_cmd
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


@register_pymol_cmd
def set_hydro_color_v2(
        selection: str = '(all)',
        palette: str = 'blue_white_red',
        scale_name: str = 'Ja',
        minimum: float = None,
        maximum: float = None,
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

    if with_sasa:
        sasa_buffer = get_sasa_by_res(selection)
    
    selection = residue_with_CA(selection)

    def set_hydro(resn, resi, chain, scale_name='Ja'):
        nonlocal minimum, maximum
        v = get_hydro(resn, scale_name=scale_name)
        residue = residue_format(resn, resi, chain, selection)

        if with_sasa:
            sasa = sasa_buffer[residue]
            v *= sasa

        if auto_bound:
            minimum = min(minimum, v)
            maximum = max(maximum, v)

        return v

    cmd.alter(
        selection,
        f'b = set_hydro(resn, resi, chain, scale_name="{scale_name}")',
        space={'set_hydro': set_hydro})
    cmd.spectrum(
        'b',
        palette=palette,
        selection=selection,
        minimum=minimum,
        maximum=maximum)


@register_pymol_cmd
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


def set_hydration(selection: str = '(all)', radius: float = 2.8, sasa_threshold: float = -1.0) -> tuple:
    """
    Count number of water molecules (hydration water molecules)
    within a given radius of each residue. These
    number will be stored in b-factor of
    each residue.

    Parameters
    ----------
    selection: str
        pymol selection string.
    radius: float
    TODO
    """
    minimum = float('inf')
    maximum = -minimum
    is_sasa = sasa_threshold >= 0

    selection = residue_with_CA(selection)

    if is_sasa:
        sasa_buffer = get_sasa_by_res(selection)

    def _update(resn, resi, chain):
        residue = residue_format(resn, resi, chain, selection)
        count = cmd.count_atoms(
            f"(resn HOH) within {radius} of ({residue} and elem N+O)")
        if is_sasa:
            if sasa_buffer[residue] < sasa_threshold:
                count = 0
        minimum = min(count, minimum)
        maximum = max(count, maximum)
        return count
    
    cmd.alter(
        selection,
        f'b = _update(resn, resi, chain)',
        space={'_update': _update})

    if minimum == float('inf'):
        raise ValueError(
            f"No valid residue detected for selection {selection}")

    return minimum, maximum


@register_pymol_cmd
def set_hydration_color(
        selection='(all)',
        radius=2.8,
        minimum: int = None,
        maximum: int = None,
        palette: str = "red_white_blue",
        sasa_threshold: float = -1.0):
    """
    Annotate residue color according
    to hydration water molecular count.

    TODO
    """
    _minimum, _maximum = set_hydration(selection, radius, sasa_threshold=sasa_threshold)
    minimum = minimum if minimum is not None else _minimum
    maximum = maximum if maximum is not None else _maximum
    cmd.spectrum('b', palette, selection, minimum=minimum, maximum=maximum)

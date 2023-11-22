from pymol import cmd
from functools import lru_cache
from ..utils import register_pymol_cmd

__all__ = ['plddt_color']


@lru_cache(maxsize=None)
def _plddt_color_defined():
    cmd.set_color('plddt_very_high', (0, 83, 214))
    cmd.set_color('plddt_confident', (101, 203, 243))
    cmd.set_color('plddt_low', (255, 219, 19))
    cmd.set_color('plddt_very_low', (255, 125, 69))


@register_pymol_cmd
def plddt_color(selection: str = '(all)'):
    """
    Plot plddt color by residue.
    Assume plddt is saved as bfactors.
    The same color as Alphafold DB.
    """
    _plddt_color_defined()
    cmd.color('plddt_very_high', selection)
    cmd.color('plddt_confident', f'b < 90 and {selection}')
    cmd.color('plddt_low', f'b < 70 and {selection}')
    cmd.color('plddt_very_low', f'b < 50 and {selection}')

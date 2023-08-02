from pymol import cmd
from . import (hydro, sasa)
from .hydro import *
from .sasa import *

__all__ = []
__all__.extend(hydro.__all__)
__all__.extend(sasa.__all__)

def __init_plugin__(*args):
    for func_name in __all__:
        cmd.extend(func_name, globals()[func_name])

from pymol import cmd
from . import (ply, dots, patch)
from .ply import *
from .dots import *
from .patch import *

__all__ = []
__all__.extend(ply.__all__)
__all__.extend(dots.__all__)
__all__.extend(patch.__all__)

def __init_plugin__(*args):
    for func_name in __all__:
        cmd.extend(func_name, globals()[func_name])


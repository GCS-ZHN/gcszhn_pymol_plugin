from pymol import cmd
from .load_ply import *
from .load_dots import *
from .extract_patch import *

def __init_plugin__(*args):
    cmd.extend('loadply', load_ply)
    cmd.extend('loaddots', load_dots)
    cmd.extend('loadgiface', load_giface)
    cmd.extend('extractpatch', extract_patch)
    cmd.extend('patchseq', patch_seq)


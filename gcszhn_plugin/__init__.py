from . import surface, pdb

def __init_plugin__(*args):
    surface.__init_plugin__(*args)
    pdb.__init_plugin__(*args)
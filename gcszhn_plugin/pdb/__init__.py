from pymol import cmd
from .hydro import set_hydro_color_v2, set_hydro_color_v1, avail_hydro_scales
from .sasa import set_sasa_color

def __init_plugin__(*args):
    cmd.extend('set_hydro_color_v1', set_hydro_color_v1)
    cmd.extend('set_hydro_color_v2', set_hydro_color_v2)
    cmd.extend('avail_hydro_scales', avail_hydro_scales)
    cmd.extend('set_sasa_color', set_sasa_color)
import logging

from pymol import cmd

# register cmd
from . import pdb as _
from . import surface as _

from .utils import __reigster_pymol_cmd__


def __init_plugin__(*args):
    for k, v in __reigster_pymol_cmd__.items():
        cmd.extend(k, v)
        logging.info(f'command {k} is registered')

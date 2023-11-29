"""
Easy tools for pymol developed
by GCS-ZHN.

See more detail at https://github.com/GCS-ZHN/gcszhn_pymol_plugin
"""

import logging
import pkgutil
from pymol import cmd

from .utils import __reigster_pymol_cmd__


def __init_plugin__(*args):
    # load all subpackage and submodule
    for module_info in pkgutil.walk_packages(__path__, __package__ + '.'):
        __import__(module_info.name)

    for k, v in __reigster_pymol_cmd__.items():
        cmd.extend(k, v)
        logging.info(f'command {k} is registered')

from pymol import cmd
from contextlib import contextmanager

__reigster_pymol_cmd__ = dict()


def register_pymol_cmd(func):
    __reigster_pymol_cmd__[func.__name__] = func
    return func


@contextmanager
def local_setting(**kwargs):
    old_settings = {k: cmd.get(k) for k in kwargs}
    for k, v in kwargs.items():
        cmd.set(k, v)
    yield kwargs
    for k, v in old_settings.items():
        cmd.set(k, v)


def residue_format(resn: str, resi: str, chain: str, selection: str) -> str:
    return f'resn {resn} and resi {resi} and chain {chain} and {selection}'


def residue_with_CA(selection: str) -> str:
    suffix = 'name CA'
    if suffix not in selection:
        selection += (' ' + suffix)
    return selection

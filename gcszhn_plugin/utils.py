from pymol import cmd
from contextlib import contextmanager


@contextmanager
def local_setting(**kwargs):
    old_settings = {k: cmd.get(k) for k in kwargs}
    for k, v in kwargs.items():
        cmd.set(k, v)
    yield kwargs
    for k, v in old_settings.items():
        cmd.set(k, v)


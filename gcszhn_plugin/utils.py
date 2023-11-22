import threading
import functools

from pymol import cmd
from typing import List
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


def lock(func):
    _lock = threading.Lock()

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with _lock:
            return func(*args, **kwargs)
    
    return wrapper


def int_array_to_str(int_array: List[int], prefix="") -> str:
    """
    Convert integer array as string.
    original array will be sorted in ascending order
    and continous value will be abbreviated
    """

    if len(int_array) == 0:
        return ""

    int_array = sorted(int_array)
    
    start = int_array[0]
    end = int_array[0]

    result = []
    for i in int_array[1:]:
        if i == end + 1:
            end = i
        else:
            if start == end:
                result.append(prefix+str(start))
            else:
                result.append(f"{prefix}{start}-{end}")
            start = i
            end = i

    if start == end:
        result.append(prefix+str(start))
    else:
        result.append(f"{prefix}{start}-{end}")

    return ",".join(result)

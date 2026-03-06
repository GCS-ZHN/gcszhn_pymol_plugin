import pandas as pd
from pymol import cmd
from pathlib import Path
from ..utils import register_pymol_cmd, int_array_to_str


@register_pymol_cmd
def copy_selection(selection: str = "(all)", mode="tab"):
    """
    Copy selection to clipboard.

    Parameters
    ----------
    selection : str, optional
        Selection.
    mode : str, optional
        Copy mode. One of "tab", "rf_hotspot", "range".
        tab: copy chain, resi, resn to clipboard in TSV format.
        rf_hotspot: copy chain and resi in RFdiffusion hotspot format.
        range: copy chain and resi in numerial range format.
    """
    def _gen():
        for atom in cmd.get_model(selection).atom:
            yield atom.chain, atom.resi, atom.resn
    
    data = pd.DataFrame(_gen(), columns=['chain', 'resi', 'resn'])
    # remove duplicated atoms
    data = data.drop_duplicates(subset=['chain', 'resi', 'resn'])
    if mode == "tab":
        data.to_clipboard(index=False)
    elif mode == "rf_hotspot":
        resi = data['chain'] + data['resi']
        resi = pd.DataFrame(resi, columns=['resi'])
        resi.T.to_clipboard(index=False, sep=",", header=False)
    elif mode == "range":
        data_group = data.groupby("chain")
        resi_range = data_group.apply(lambda d: int_array_to_str(d['resi'].astype(int), d.name))
        resi_range.to_clipboard(index=False, header=False)
    else:
        raise ValueError(f"Unknown mode: {mode}")


@register_pymol_cmd
def save_by_objects(prefix: str = "", file_type: str = "pdb"):
    """
    Save objects to PDB files.

    Parameters
    ----------
    prefix : str, optional
        Prefix of PDB files.
    file_type: str, optional
        Any supported file type in pymol.
    """
    for obj in cmd.get_object_list():
        cmd.save(f"{prefix}{obj}.{file_type}", obj)


@register_pymol_cmd
def png_pse(width: int, height: int, transparent: bool = False, output: str = None) -> str:
    """
    Save PNG image with default pse path.

    Parameters
    ----------
    width : int
        Image width.
    height : int
        Image height.
    transparent : bool, optional
        Whether to save image with transparent background.
    output : str, optional
        Output path. If not specified, the image will be saved
        to the same directory as the session file with the same name
        and .png extension.
    """
    if output is None:
        session_file = cmd.get("session_file")
        if not session_file.strip():
            raise ValueError("No session file found. Please specify output path.")
        output = str(Path(session_file).with_suffix(".png"))
    
    width = int(width)
    height = int(height)
    if transparent:
        cmd.set("ray_opaque_background", 'off')
    cmd.ray(width, height)
    cmd.png(output)
    print(f'png_pse: saved image to {output} (width={width}, height={height}, transparent={transparent})')
    return output

# GCSZHN Pymol plugin

download and packaged "gcszhn_plugin" as a ZIP file, and install it as pymol plugin.

## provided pymol command

- copy_selection

Copy selection to clipboard.

- save_by_objects

Save different objects as different files.

- plddt_color

Plot plddt color by residue. The same color as Alphafold DB.

- set_hydro_color_v2

Annotate hydrophobicity color of each residue. this function implmented by `pymol.cmd.alter` and `pymol.cmd.spectrum`.

- set_hydro_color_v1

Annotate hydrophobicity color of each residue.

- set_hydration_color

Annotate residue color according to hydration water molecular count.

- get_sasa

Calculate solvent accessible surface area (SASA) for specific selection.

- set_sasa_color

Annotated color by sasa.

- load_dots

Load 3D object (*.ply) as dots.

- load_ply

Load 3D object (*.ply) as mesh.

- load_ply_with_patch

Load 3D object (*.ply) as mesh and annotate specific patch by provided patch list.

- load_giface
Load 3D object (*.ply) as giface.

- extract_patch

Extract patch residues from model according to the patch mesh object.

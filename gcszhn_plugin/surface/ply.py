from pymol import cmd
from colour import Color
from itertools import combinations
from pymol.cgo import *
from .mesh_utils import Mesh
import numpy as np

__all__ = ['load_ply', 'load_giface']

colorDict = {'sky': [0.0, 0.76, 1.0 ],
        'sea': [0.0, 0.90, 0.5 ],
        'yellowtint': [0.88, 0.97, 0.02 ],
        'hotpink': [0.90, 0.40, 0.70 ],
        'greentint': [0.50, 0.90, 0.40 ],
        'blue': [0.0, 0.0, 1.0 ],
        'green': [0.0, 1.0, 0.0 ],
        'yellow': [1.0, 1.0, 0.0 ],
        'orange': [1.0, 0.5, 0.0],
        'red': [1.0, 0.0, 0.0],
        'black': [0.0, 0.0, 0.0],
        'white': [1.0, 1.0, 1.0],
        'gray': [0.9, 0.9, 0.9] }

# Create a gradient color from color 1 to whitish, to color 2. val goes from 0 (color1) to 1 (color2).
def color_gradient(vals, color1, color2):
    c1 = Color("white")
    c2 = Color("orange")
    ix = np.floor(vals*100).astype(int)
    crange = list(c1.range_to(c2, 100))
    mycolor = []
    print(crange[0].get_rgb())
    for x in ix: 
        myc = crange[x].get_rgb()
        mycolor.append([myc[0], myc[1], myc[2]]) 
    return mycolor

def iface_color(iface):
    # max value is 1, min values is 0
    hp = iface.copy()
    hp = hp*2 - 1
    mycolor = charge_color(-hp)
    return mycolor

# Returns the color of each vertex according to the charge. 
# The most purple colors are the most hydrophilic values, and the most 
# white colors are the most positive colors.
def hphob_color(hphob):
    # max value is 4.5, min values is -4.5
    hp = hphob.copy()
    # normalize
    hp = hp + 4.5 
    hp = hp/9.0
    #mycolor = [ [COLOR, 1.0, hp[i], 1.0]  for i in range(len(hp)) ]
    mycolor = [ [1.0, 1.0-hp[i], 1.0]  for i in range(len(hp)) ]
    return mycolor

# Returns the color of each vertex according to the charge. 
# The most red colors are the most negative values, and the most 
# blue colors are the most positive colors.
def charge_color(charges):
    # Assume a std deviation equal for all proteins.... 
    max_val = 1.0
    min_val = -1.0

    norm_charges = charges
    blue_charges = np.array(norm_charges)
    red_charges = np.array(norm_charges)
    blue_charges[blue_charges < 0] = 0
    red_charges[red_charges > 0] = 0
    red_charges = abs(red_charges) 
    red_charges[red_charges>max_val] = max_val
    blue_charges[blue_charges< min_val] = min_val
    red_charges = red_charges/max_val
    blue_charges = blue_charges/max_val
    #red_charges[red_charges>1.0] = 1.0
    #blue_charges[blue_charges>1.0] = 1.0
    green_color  = np.array([0.0]*len(charges))
    mycolor = [ [0.9999-blue_charges[i], 0.9999-(blue_charges[i]+red_charges[i]), \
                    0.9999-red_charges[i]]  for i in range(len(charges)) ]
    for i in range(len(mycolor)):
        for k in range(3):
            if mycolor[i][k] < 0:
                mycolor[i][k] = 0

    return mycolor


def si_color(si):
    return charge_color(si)


def ddc_color(ddc):
    return charge_color(ddc)

def hbond_color(hbond):
    return charge_color(hbond)


def add_triangle_faces(faces, vertices, colors, normals):
    obj = []
    for triangle in faces:
        obj.extend([BEGIN, TRIANGLES])
        for i in range(3):
            i = int(i)
            obj.append(COLOR)
            obj.extend(colors[triangle[i]])
            obj.append(NORMAL)
            obj.extend(normals[triangle[i]])
            obj.append(VERTEX)
            obj.extend(vertices[triangle[i]])
        obj.extend([END])
    return obj


def load_ply_with_patch(ply_file, patch_list_file, *patch_id_colors, name = None, backgroud_color = 'gray'):
    if not name:
        name = os.path.basename(ply_file).split('.')[0]
    mesh = Mesh.create_mesh()
    mesh.load_mesh(ply_file)
    patch_list = np.load(patch_list_file, allow_pickle=True)
    patch_dict = {patch[0]: patch for patch in patch_list}
    vertices = mesh.vertices
    faces = mesh.faces
    colors = [colorDict[backgroud_color]]*len(vertices)
    for patch_id, color in patch_id_colors:
        patch = patch_dict[patch_id]
        for i in patch:
            if isinstance(color, str):
                colors[i] = colorDict[color]
            elif isinstance(color, list):
                colors[i] = color
            else:
                raise ValueError("Color must be a string or a list of 3 values")
    nx = mesh.get_attribute('vertex_nx')
    ny = mesh.get_attribute('vertex_ny')
    nz = mesh.get_attribute('vertex_nz')
    normals = np.vstack([nx, ny, nz]).T

    obj = add_triangle_faces(faces, vertices, colors, normals)
    cmd.load_cgo(obj, name)



def load_ply(filename, group_name = None, vertex_size=0.2, enable_properties = None):
    mesh = Mesh.create_mesh()
    mesh.load_mesh(filename)
    if not group_name:
        group_name = os.path.basename(filename).split('.')[0]
    group_members = []
    with_normal = False
    verts = mesh.vertices
    faces = mesh.faces

    try:
        if enable_properties is None or enable_properties[0] == 'vertex_charge':
            color_array = charge_color(mesh.get_attribute("vertex_charge"))
        elif enable_properties[0] == 'vertex_hphob':
            color_array = hphob_color(mesh.get_attribute("vertex_hphob"))
        elif enable_properties[0] == 'vertex_si':
            color_array = si_color(mesh.get_attribute("vertex_si"))
        elif enable_properties[0] == 'vertex_ddc':
            color_array = ddc_color(mesh.get_attribute('vertex_ddc') * 1.4285)
        elif enable_properties[0] == 'vertex_iface':
            color_array = iface_color( mesh.get_attribute("vertex_iface"))
        elif enable_properties[0] == 'vertex_hbond':
            color_array = hbond_color(mesh.get_attribute("vertex_hbond"))
        else:
            color_array = [colorDict['green']]*len(verts)
    except:
        color_array = [colorDict['green']]*len(verts)

    if 'vertex_nx' in mesh.get_attribute_names():
        nx = mesh.get_attribute('vertex_nx')
        ny = mesh.get_attribute('vertex_ny')
        nz = mesh.get_attribute('vertex_nz')
        normals = np.vstack([nx, ny, nz]).T
        with_normal = True

    # Draw vertices 
    obj = []
    for v_ix in range(len(verts)):
        vert = verts[v_ix]
        colorToAdd = color_array[v_ix]
        # Vertices
        obj.append(COLOR)
        obj.extend(colorToAdd)
        obj.extend([SPHERE, vert[0], vert[1], vert[2], vertex_size])

    group_members.append(group_name + "_vertices")
    cmd.load_cgo(obj, group_members[-1], 1.0)


    if with_normal:
        # Draw surface charges.
        if 'vertex_charge' in mesh.get_attribute_names() and (enable_properties is None or 'vertex_charge' in enable_properties): 
            color_array_surf = charge_color(mesh.get_attribute("vertex_charge"))
            obj = add_triangle_faces(faces, verts, color_array_surf, normals)
            group_members.append(group_name + "_vertex_charge")
            cmd.load_cgo(obj, group_members[-1], 1.0)

        # Draw hydrophobicity
        if 'vertex_hphob' in mesh.get_attribute_names() and (enable_properties is None or 'vertex_hphob' in enable_properties): 
            hphob = mesh.get_attribute('vertex_hphob')
            color_array_surf = hphob_color(hphob)
            obj = add_triangle_faces(faces, verts, color_array_surf, normals)
            group_members.append(group_name + "_vertex_hphob")
            cmd.load_cgo(obj, group_members[-1], 1.0)

        # Draw shape index
        if 'vertex_si' in mesh.get_attribute_names() and (enable_properties is None or 'vertex_si' in enable_properties): 
            color_array_surf = si_color(mesh.get_attribute('vertex_si'))
            obj = add_triangle_faces(faces, verts, color_array_surf, normals)
            group_members.append(group_name + "_vertex_si")
            cmd.load_cgo(obj, group_members[-1], 1.0)

        # Draw ddc
        if 'vertex_ddc' in mesh.get_attribute_names() and (enable_properties is None or 'vertex_ddc' in enable_properties): 
            # Scale to -1.0->1.0
            color_array_surf = ddc_color(mesh.get_attribute('vertex_ddc') * 1.4285)
            obj = add_triangle_faces(faces, verts, color_array_surf, normals)
            group_members.append(group_name + "_vertex_ddc")
            cmd.load_cgo(obj,group_members[-1], 1.0)

        # Draw iface
        if 'vertex_iface' in mesh.get_attribute_names() and (enable_properties is None or 'vertex_iface' in enable_properties): 
            color_array_surf = iface_color(mesh.get_attribute('vertex_iface'))
            obj = add_triangle_faces(faces, verts, color_array_surf, normals)
            group_members.append(group_name + "_vertex_iface")
            cmd.load_cgo(obj, group_members[-1], 1.0)

        # Draw hbond
        if 'vertex_hbond' in mesh.get_attribute_names() and (enable_properties is None or 'vertex_hbond' in enable_properties): 
            color_array_surf = charge_color(mesh.get_attribute('vertex_hbond'))
            obj = add_triangle_faces(faces, verts, color_array_surf, normals)
            group_members.append(group_name + "_vertex_hbond")
            cmd.load_cgo(obj, group_members[-1], 1.0)

        # Draw normals
        obj = []
        if enable_properties is None or 'normal' in enable_properties:
            for v_ix in range(len(verts)):
                colorToAdd = colorDict['white']
                vert1 = verts[v_ix]
                vert2 = [verts[v_ix][0]+nx[v_ix],\
                        verts[v_ix][1]+ny[v_ix],\
                        verts[v_ix][2]+nz[v_ix]]
                obj.extend([LINEWIDTH, 2.0])
                obj.extend([BEGIN, LINES])
                obj.append(COLOR)
                obj.extend(colorToAdd)
                obj.extend([VERTEX, (vert1[0]), (vert1[1]), (vert1[2])])
                obj.extend([VERTEX, (vert2[0]), (vert2[1]), (vert2[2])])
                obj.append(END)
            group_members.append(group_name + "_normal")
            cmd.load_cgo(obj, group_members[-1], 1.0)


    # Draw triangles (faces)
    if enable_properties is None or 'mesh' in enable_properties:
        obj = []
        for tri in faces: 
            colorToAdd = colorDict['gray']
            for line in combinations(tri, 2): 
                obj.extend([BEGIN, LINES])
                obj.append(COLOR)
                obj.extend(colorToAdd)
                for vert in line: 
                    obj.append(VERTEX)
                    obj.extend(verts[vert])
                obj.append(END)
        group_members.append(group_name + "_mesh") 
        cmd.load_cgo(obj, group_members[-1], 1.0)

    cmd.group(group_name, " ".join(group_members))

# Load the sillouete of an iface.
def load_giface(filename, color="white", name='giface', dotSize=0.2, lineSize = 1.0):
    mesh = Mesh.create_mesh()
    mesh.load_mesh(filename)
    if 'vertex_iface' not in mesh.get_attribute_names():
        return
    iface = mesh.get_attribute('vertex_iface')
    # Color an edge only if:
        # iface > 0 for its two edges
        # iface is zero for at least one of its edges.
    # Go through each face. 
    faces = mesh.faces
    verts = mesh.vertices
    obj = []
    visited = set()
    colorToAdd = colorDict['green']
    obj.extend([BEGIN, LINES])
    obj.extend([LINEWIDTH, 5.0])
    obj.extend(colorToAdd)
    for tri in faces: 
        pairs = [[tri[0],tri[1], tri[2]], [tri[0],tri[2], tri[1]], [tri[1],tri[2], tri[0]]]
        for pair in pairs: 
            if iface[pair[0]] > 0 and iface[pair[1]] > 0 and iface[pair[2]] == 0:
                vert1 = verts[pair[0]]
                vert2 = verts[pair[1]]

                obj.extend([VERTEX, (vert1[0]), (vert1[1]), (vert1[2])])
                obj.extend([VERTEX, (vert2[0]), (vert2[1]), (vert2[2])])
    obj.append(END)
    name = "giface_"+filename 
    cmd.load_cgo(obj,name, 1.0)
    colorToAdd = colorDict['green']

    obj = []
    obj.extend(colorToAdd)
    for tri in faces: 
        pairs = [[tri[0],tri[1], tri[2]], [tri[0],tri[2], tri[1]], [tri[1],tri[2], tri[0]]]
        for pair in pairs: 
            if iface[pair[0]] > 0 and iface[pair[1]] > 0 and iface[pair[2]] == 0:
                vert1 = verts[pair[0]]
                vert2 = verts[pair[1]]

                obj.extend([SPHERE, (vert1[0]), (vert1[1]), (vert1[2]), 0.4])
                obj.extend([SPHERE, (vert2[0]), (vert2[1]), (vert2[2]), 0.4])
    #obj.append(END)
    name = "giface_verts_"+filename 
    cmd.load_cgo(obj,name, 1.0)
 

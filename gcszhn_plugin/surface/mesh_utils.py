from typing import Iterable
import numpy as np
from abc import ABCMeta, abstractmethod


class Mesh(metaclass=ABCMeta):
    """Abstract base class for mesh objects."""

    def __init__(self):
        self.vertices = None
        self.faces = None

    @abstractmethod
    def load_mesh(self, filename: str):
        """Load mesh from file."""
        raise NotImplementedError
    
    @abstractmethod
    def get_attribute_names(self) -> Iterable[str]:
        """Return a list of attribute names."""
        raise NotImplementedError
    
    @abstractmethod
    def get_attribute(self, attribute_name: str) -> np.ndarray:
        """Return a copy of the attribute."""
        raise NotImplementedError
    
    @staticmethod
    def create_mesh():
        if MESHIO_AVAILABLE:
            return MeshioAdaptor()
        elif PYMESH_AVAILABLE:
            return PyMeshAdaptor()
        else:
            return SimpleMesh()


try:
    import meshio
    MESHIO_AVAILABLE = True

    class MeshioAdaptor(Mesh):
        """Meshio adaptor for triangle mesh."""

        def __init__(self):
            super().__init__()
            self.mesh = None

        def load_mesh(self, filename: str):
            self.mesh = meshio.read(filename)
            self.vertices = self.mesh.points
            self.faces = self.mesh.cells_dict['triangle']

        def get_attribute_names(self)-> Iterable[str]:
            for key in self.mesh.point_data.keys():
                yield 'vertex_' + key
            
            if 'triangle' in self.mesh.cell_data:
                for key in self.mesh.cell_data['triangle'].keys():
                    yield 'face_' + key

        def get_attribute(self, attribute_name: str)-> np.ndarray:
            if attribute_name.startswith('vertex_'):
                attribute_name = attribute_name[7:]
                return np.copy(self.mesh.point_data[attribute_name])
            elif attribute_name.startswith('face_'):
                attribute_name = attribute_name[5:]
                return np.copy(self.mesh.cell_data['triangle'][attribute_name])

except ImportError:
    MESHIO_AVAILABLE = False

try:
    import pymesh
    PYMESH_AVAILABLE = True
    class PyMeshAdaptor(Mesh):
        """PyMesh adaptor for triangle mesh."""

        def __init__(self):
            super().__init__()
            self.mesh = None

        def load_mesh(self, filename: str):
            self.mesh = pymesh.load_mesh(filename)
            self.vertices = self.mesh.vertices
            self.faces = self.mesh.faces

        def get_attribute_names(self)-> Iterable[str]:
            return self.mesh.get_attribute_names()

        def get_attribute(self, attribute_name: str)-> np.ndarray:
            return np.copy(self.mesh.get_attribute(attribute_name))
except ImportError:
    PYMESH_AVAILABLE = False


class SimpleMesh(Mesh):
    """
    Simple mesh class to load ply files."""
    def __init__(self):
        super().__init__()
        self.vertices = []
        self.faces = []

    def load_mesh(self, filename: str):
        lines = open(filename, 'r').readlines()
        # Read header
        self.attribute_names = []
        self.num_verts = 0
        line_ix = 0
        while 'end_header' not in lines[line_ix]: 
            line = lines[line_ix]
            if line.startswith('element vertex'): 
                self.num_verts = int(line.split(' ')[2])
            if line.startswith('property float') or line.startswith('property double'):
                self.attribute_names.append('vertex_'+line.split(' ')[2].rstrip())
            if line.startswith('element face'):
                self.num_faces= int(line.split(' ')[2])
            line_ix += 1
        line_ix += 1
        header_lines = line_ix
        self.attributes = {}
        for at in self.attribute_names:
            self.attributes[at] = []
        self.vertices = []
        self.normals = []
        self.faces = []
        # Read vertex attributes.
        for i in range(header_lines, self.num_verts+header_lines):
            cur_line = lines[i].split(' ')
            vert_att = [float(x) for x in cur_line]
            # Organize by attributes
            for jj, att in enumerate(vert_att): 
                self.attributes[self.attribute_names[jj]].append(att)
            line_ix += 1
        # Set up vertices
        for jj in range(len(self.attributes['vertex_x'])):
            self.vertices = np.vstack([self.attributes['vertex_x'],\
                                    self.attributes['vertex_y'],\
                                    self.attributes['vertex_z']]).T
        # Read faces.
        face_line_start = line_ix
        for i in range(face_line_start, face_line_start+self.num_faces):
            try:
                fields = lines[i].split(' ')
            except:
                print('Error reading line %d' % i)
                raise
            face = [int(x) for x in fields[1:]]
            self.faces.append(face)
        self.faces = np.array(self.faces)
        self.vertices = np.array(self.vertices)
        # Convert to numpy array all attributes.
        for key in self.attributes.keys():
            self.attributes[key] = np.array(self.attributes[key])

    def get_attribute_names(self) -> Iterable[str]:
        for key in self.attribute_names:
            yield key

    def get_attribute(self, attribute_name: str) -> np.ndarray:
        return np.copy(self.attributes[attribute_name])

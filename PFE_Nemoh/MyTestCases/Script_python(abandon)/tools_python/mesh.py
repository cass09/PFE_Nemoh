from math import pi
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
import shutil


class SurF(object) :

    def __init__(self, coord) :
        """
        coord (3D numpy array) : (Number of panels,
                                  Number of nodes for each panel,
                                  Number of coordinates for each node == 3)
        """
        self.coord = coord
        
    def info(self):
        """
        """
        Npanels = self.coord.shape[0]
        nodes = list()
        conectivity = np.zeros((Npanels, 4), dtype = int)
        for i0, nodeori in enumerate(self.coord.reshape((-1,3))) :
            nodes.append(nodeori)
            for i1, node in enumerate(nodes) :
                if all(abs(nodeori-node) < 1e-6) :
                    conectivity[i0/4, int((i0/4.-i0/4)*4)] = i1+1
                    if i1+1 < len(nodes) :
                        nodes.pop()
                        break
        nodes = np.array(nodes, float)
        return nodes, conectivity-1

    

    def translation(self, disp) :
        """
        Translation is applied to coord according to
        a the displacement given by disp.

        disp (1D numpy array) : [translation in x,
                                 translation in y,
                                 translation in z]
        """
        self.coord += disp

    def rotation(self, rot) :
        """
        Rotation is applied to coord according to
        rot. The rotation is performend with respect to the (0,0,0)
        point in self.coord.

        rot (1D numpy array) : [rotation in x,
                                rotation in y,
                                rotation in z] (radians)
        """
        nodes = self.coord.reshape((-1,3))
        R = [np.array([[np.cos(angl),-np.sin(angl)],
                       [np.sin(angl),np.cos(angl)]],
                        dtype = float) for angl in rot]
        for ind0 in range(len(nodes)) :
            for ind1, ax in enumerate(((1,2),(0,2),(0,1))) :
                nodes[ind0, ax] = np.dot(R[ind1],nodes[ind0, ax])
        self.coord = nodes.reshape((-1,4,3))
  

    def dat(self, fn, sym = 0) :
        """
        Generate a .dat file for the surface

        fn (string) : direction and name of the .dat file that will be
                 generated. e.g. ".\\DesiredFolder\\DesiredName.dat"
        refinement (bool) : if True (default), the format will follow that
                            of Nemoh mesh file (the one given in Nemoh.cal).
                            If False, the format will follow that prior
                            refinement for mesh generation using Mesh.exe.
        sym (int): if refined then sym is 1 if a symmetry about the (xOz)
                   plane is used. 0 otherwise
        """
        fn = fn.split('.dat')[0]
        fn += '.dat'
        Npanels = self.coord.shape[0]
        nodes = list()
        conectivity = np.zeros((Npanels, 4), dtype = int)
        for i0, nodeori in enumerate(self.coord.reshape((-1,3))) :
            nodes.append(nodeori)
            for i1, node in enumerate(nodes) :
                if all(nodeori == node) :
                    conectivity[int(i0/4), int((i0/4.-i0//4)*4)] = i1+1
                    if i1+1 < len(nodes) :
                        nodes.pop()
                        break
        with open(fn , 'w') as fid:
            fid.write('2 {:}\n'.format(sym))
            for i0, node in enumerate(nodes) :
                fid.write('{:} '.format(i0+1))
                fid.write('{:} {:} {:}\n'.format(*node))
            fid.write('0 0. 0. 0.\n')
            for panel in conectivity :
                fid.write('{:} {:} {:} {:}\n'.format(*panel))
            fid.write('0 0. 0. 0.\n')
        
        self.dirDAT = fn
        self.Npanels = Npanels
        self.Nnodes = len(nodes)

    def save_mesh(self, original_path, new_path, change=False):
        """
        Save file mesh.dat.
        If file changed, file rewritten.
        Else just copy past.
        """
        if change:
            self.dat(new_path) 
        else:
            shutil.copy(original_path, new_path)
            self.dirDAT = new_path  
            self.Npanels = self.coord.shape[0] 

            unique_nodes = set()
            for panel in self.coord:
                for node in panel:
                    unique_nodes.add(tuple(node))
            self.Nnodes = len(unique_nodes)


def readDAT(fn, translation = np.zeros(3), rotation = np.zeros(3)) :
    """
    Reads .dat mesh files and returns the x,y,z coordinates of the nodes
    of each panel through an instance of SurF. If translation or rotation
    are non zero arrays, a translated and/or rotated SurF
    instance will be generated instead.

    fn (string) : direction and name of the .dat file that will be
                 generated. e.g. ".\\DesiredFolder\\DesiredName.dat"
    translation (1D numpy array) : [translation in x,
                                    translation in y,
                                    translation in z]
    rotation (1D numpy array) : [rotation in x,
                                 rotation in y,
                                 rotation in z]
                                 with respect to the (0,0,0)
                                 point given in the mesh file. In this
                                 regard, it is worth noticing that
                                 the translation is performed before
                                 the rotation.
    """
    ## read dat
    with open(fn,'r') as fIlE :
        nodes = list()
        conectivity = list()
        line = fIlE.readline().split('!')[0].split()
        if int(line[0]) == 2 :
            sym = int(line[1])
        else :
            Nn = int(line[0])
            sym = -1
        if sym > -1 :
            trigger = False
            for fl in fIlE.readlines() :
                line = np.array(fl.split('!')[0].split(), dtype = float)
                if not trigger :
                    nodes.append(line)
                elif trigger :
                    conectivity.append(line)
                if all(line == 0) :
                    if not trigger:
                        nodes.pop()
                    elif trigger :
                        conectivity.pop()
                    trigger = True
            nodes = np.array(nodes, dtype = float).reshape((-1,4))[:,1:]
        else :
            raise ValueError("Format file incorrect - Nemoh mesh.dat expected")
        conectivity = np.array(conectivity, dtype = int).reshape((-1,4))-1
        coord = np.zeros((len(conectivity),4,3), dtype = float)
        for panel , i0 in enumerate(conectivity):
            for node in range(4) :
                coord[panel,node] = nodes[i0[node]]
                if any(translation != 0) :
                    coord[panel,node] += translation
                if any(rotation != 0) :
                    for d,ax in enumerate(((1,2),(0,2),(0,1))) :
                        R = np.array([[np.cos(rotation[d]),-np.sin(rotation[d])],
                                      [np.sin(rotation[d]),np.cos(rotation[d])]], dtype = float)
                        coord[panel,node,ax] = np.dot(R,[panel,node,ax])
    ## instantiate the surface using coord as (Npanels, Nnodes, Ncoord)
    surface = SurF(coord)
    return surface

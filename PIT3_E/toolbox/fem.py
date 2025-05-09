# -*- coding: utf-8 -*-
"""
"""

import numpy as np
import matplotlib.pyplot as plt

class fem(object):
    """
    """
    
    def __init__(self, nodes, conectivity):
        """
        Only 2D triangle elements consisting of 3 nodes, piecewise solution representation and
        linear polynomial inside-element interpolation shape functions.
        
        nodes.shape = (Nnodes, 2) 1st column x-coord and 2nd y-coord
        conectivity.shape = (Nelements, 3) 1st column id first node and so on and so forth
        """
        
        self.nodes = nodes
        self.conectivity = conectivity
        self.nodesIN = np.arange(len(nodes))
        self.phi = np.zeros(len(nodes), complex) # trivial solution
    
    def field2D(self, A, B, f):
        """
        2D divergence and gradient operators
        
        div(A*grad(phi)) + B*phi = f
        
        A.shape = (Nelements, ) A is constant within an element
        B.shape = (Nelements, ) B is constant within an element
        f.shape = (Nelements, ) f is constant within an element
        """
        
        Nnodes = len(self.nodes)
        Bsys = np.zeros(Nnodes, complex)
        Msys = np.zeros((Nnodes, Nnodes), complex)
        for el, cnctel in enumerate(self.conectivity):
            xl, xm, xn = self.nodes[cnctel, 0]
            yl, ym, yn = self.nodes[cnctel, 1]
            dJ = (xl-xn)*(ym-yn)-(xm-xn)*(yl-yn) # determinant of jacobian matrix
            T = 1./dJ*np.array([[ym-yn, -(xm-xn)], 
                                [-(yl-yn), xl-xn]], float) # x_local = T.dot(x_global)
            gradN = np.array([[T[0,0],           T[0,1]],
                              [T[1,0],           T[1,1]],
                              [-(T[0,0]+T[1,0]), -(T[0,1]+T[1,1])]], float)
            
            msys = -1./24.*np.ones((3, 3), float)*B[el]*np.abs(dJ)
            msys[np.arange(3), np.arange(3)] *= 2.
            msys += 1./2.*gradN.dot(gradN.T)*A[el]*np.abs(dJ)
              
            eq, unk = np.meshgrid(cnctel, cnctel, indexing='ij')
            Msys[eq, unk] += msys

            Bsys[eq] += -1./6.*np.ones(3, float)*f[el]*np.abs(dJ)
        
        self.Msys = Msys
        self.Bsys = Bsys
        
    def Dirichlet(self, nodesDirichlet, phiDirichlet):
        """
        nodesDirichlet.shape = (Nnodes in Dirichlet boundary, )
        phiDirichlet.shape = (Nnodes in Dirichlet boundary,) the corresponding prescribed phi
        """
        
        # inclusion of the Dirichlet boundary conditions
        for nde, phi in zip(nodesDirichlet, phiDirichlet):
            self.Bsys = self.Bsys[self.nodesIN!=nde] # remove equations
            self.Msys = self.Msys[self.nodesIN!=nde] # remove equations
            self.Bsys += -self.Msys[:, self.nodesIN==nde].reshape(-1)*phi # move known phi to the right-hand side
            self.Msys = self.Msys[:, self.nodesIN!=nde] # remove unknowns
            self.nodesIN = self.nodesIN[self.nodesIN!=nde] # remove node
        
        self.nodesDirichlet = nodesDirichlet
        self.phi[nodesDirichlet] = phiDirichlet
        
    def Neumann(self, ANeumann, conectNeumann, dphidnNeumann):
        """
        div(A*grad(phi)) + B*phi = f --> ANeumann
        d(phi)/d(n) = dphidnNeumann for all (x,y) in Neumann boundary
        
        ANeumann.shape = (Nelements in Neumann boundary,)
        conectNeumann.shape = (Nelements in Neumann boundary, )
        dphidnNeumann.shape = (Nelements in Neumann boundary,) the corresponding prescribed d(phi)/d(n),
        considered constant along the corresponding boundary element
        """
        
        # inclusion of the Neumann boundary conditions
        for el, cnctel in enumerate(conectNeumann):
            dl = np.sqrt(((self.nodes[cnctel[0]]-self.nodes[cnctel[1]])**2).sum())
            cond = (self.nodesIN==cnctel[0])+(self.nodesIN==cnctel[1])
            self.Bsys[cond] += ANeumann[el]*dphidnNeumann[el]*1./2.*dl

    def Robin(self, ARobin, conectRobin, C, g):
        """
        div(A*grad(phi)) + B*phi = f --> ARobin
        d(phi)/d(n) + C*phi = g for all (x,y) in Robin boundary
        
        ARobin.shape = (Nelements in Robin boundary,)
        conectRobin.shape = (Nelements in Robin boundary, )
        C.shape = (Nelements in Robin boundary,) considered constant along the corresponding boundary element
        g.shape = (Nelements in Robin boundary,) considered constant along the corresponding boundary element
        """
        
        # inclusion of the Neumann boundary conditions
        for el, cnctel in enumerate(conectRobin):
            dl = np.sqrt(((self.nodes[cnctel[0]]-self.nodes[cnctel[1]])**2).sum())
            cond = (self.nodesIN==cnctel[0])+(self.nodesIN==cnctel[1])
            self.Bsys[cond] += ARobin[el]*g[el]*1./2.*dl
            
            ind = np.arange(len(self.nodesIN))[cond]
            eq, unk = np.meshgrid(ind, ind, indexing='ij')
            self.Msys[eq, unk] += ARobin[el]*C[el]*np.array([[1./3., 1./6.], [1./6., 1./3.]])*dl
            
        
    def solve(self):
        """
        """
     
        self.phi[self.nodesIN] = np.linalg.solve(self.Msys, self.Bsys)
        
    def show(self):
        """
        """
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        tcl = ax.tripcolor(self.nodes[:, 0], self.nodes[:, 1], np.abs(self.phi))
        ax.plot(self.nodes[self.nodesDirichlet, 0], self.nodes[self.nodesDirichlet, 1], 'sk-')
        Dx, Dy = self.nodes[:, 0].max()-self.nodes[:, 0].min(), self.nodes[:, 1].max()-self.nodes[:, 1].min()
        ax.set_xlim(self.nodes[:, 0].min()-Dx*.002, self.nodes[:, 0].max()+Dx*.002)
        ax.set_ylim(self.nodes[:, 1].min()-Dy*.002, self.nodes[:, 1].max()+Dy*.002)
        ax.set_xlabel('x', fontsize=20)
        ax.set_ylabel('y', fontsize=20)
        cb = plt.colorbar(tcl)
        cb.set_label('abs(phi)', fontsize=20)
        
        return fig

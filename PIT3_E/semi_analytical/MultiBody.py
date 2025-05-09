"""
"""

# Imports
from math import pi
import numpy as np
from numpy.linalg import solve
from scipy.special import jv, yv, kv
from toolbox.CalaixSastre import len2
from glob import glob
import os
# Settings


class MultiBody(object):

    def __init__(self, directions, iBody, cylamplitude= True):
        """
        """
        # cylamplitude = True ==> will save cylindrical amplitude coefficients
        self.cylamplitude = cylamplitude
        # wave data from iBody
        self.depth = iBody.depth
        self.period = iBody.period
        self.wnumber = iBody.wnumber
        self.wnumber_E = iBody.wnumber_E
        # wave directions
        self.dir = directions
        # Body instance
        self.Body = iBody

    def reduce_num_periods(self, per_red):
        """
        """
        cond = np.array([np.arange(len(self.period))[self.period == per][0] for per in per_red], int)
        self.period = self.period[cond]
        self.wnumber = self.wnumber[cond]
        self.wnumber_E = self.wnumber_E[cond]
        self.Body.Fex = self.Body.Fex[cond]
        self.Body.Madd = self.Body.Madd[cond]
        self.Body.Crad = self.Body.Crad[cond]
        self.Body.D = self.Body.D[cond]
        self.Body.G = self.Body.G[cond]
        self.Body.AR = self.Body.AR[cond]
        self.Body.period = self.Body.period[cond]
        self.Body.wnumber = self.Body.wnumber[cond]
        self.Body.wnumber_E = self.Body.wnumber_E[cond]
        self.Body.truncorder = self.Body.truncorder[cond]

    def Scattering(self, coord, Evanescent):
        """
        """
        Nb = len2(coord)
        Nm = self.Body.order.max()
        dim = 2*Nm+1
        k0 = self.wnumber
        kl = self.wnumber_E
        direction = self.dir
        #
        T = self.Transformation(coord, Nm, Evanescent)
        T, D_array, G_array, M_int = self.Interaction(T, Nm, 'S', Evanescent)
        # Get the ambient planar wave, AP, for the scattering problem for the entire array
        matk0,matdir,X,mode = np.meshgrid(k0, direction, coord[:,0], range(-Nm,Nm+1),indexing='ij')
        Y = np.meshgrid(k0, direction, coord[:,1], range(-Nm,Nm+1), indexing='ij')[2]
        # Calculate AP using (k0,dir,bodycoord,mode)
        if self.Body.convention == 'N' :
            if Evanescent :
                Ntot=(2*Nm+1)*(len2(kl[0])+1)
                AP = np.zeros((len2(k0),len2(direction),len2(coord), Ntot), dtype=complex)
                for w in range(len2(k0)) : 
                    for beta in range(len2(direction)):
                        for body in range(len2(coord)):
                            for m in range(dim):
                                AP[w, beta, body, m] = np.exp(1j*mode[w, beta, body, m]*(pi/2-matdir[w, beta, body, m]))
                                AP[w, beta, body, m] *= np.exp(1j*matk0[w, beta, body, m]*(X[w, beta, body, m]*np.cos(matdir[w, beta, body, m])+Y[w, beta, body, m]*np.sin(matdir[w, beta, body, m])))
                                for l in range(len2(kl[0])):
                                    AP[w, beta, body, (2*Nm+1)*(l+1)+m] = np.exp(1j*mode[w, beta, body, m]*(pi/2-matdir[w, beta, body, m]))
                                    AP[w, beta, body, (2*Nm+1)*(l+1)+m] *= np.exp(1j*kl[w][l]*(X[w, beta, body, m]*np.cos(matdir[w, beta, body, m])+Y[w, beta, body, m]*np.sin(matdir[w, beta, body, m])))
                                    # AP[w, beta, body, (2*Nm+1)*(l+1)+m] = 0
            else : 
                AP = np.exp(1j*mode*(pi/2-matdir))
                AP *= np.exp(1j*matk0*(X*np.cos(matdir)+Y*np.sin(matdir)))
        else :
            AP = np.exp(-1j*mode*(pi/2+matdir))
            AP *= np.exp(-1j*matk0*(X*np.cos(matdir)+Y*np.sin(matdir)))
        if Evanescent : 
            AP = np.reshape(AP,(len2(k0),len2(direction),(2*Nm+1)*Nb*(len2(kl[0])+1)))
        else :     
            AP = np.reshape(AP,(len2(k0),len2(direction),(2*Nm+1)*Nb))
        Fex = np.zeros((len2(k0),len2(direction),len2(G_array[0][0])), dtype=complex)
        # Save amplitude coefficients
        if self.cylamplitude:
            self.aS = []
            self.AP = []
        for i in range(len2(k0)):
            Nmi = (len2(T[i])//Nb-1)//2
            dimi = 2*Nmi+1
            col = (np.array([range(dimi)]*Nb).T+np.array(range(0,dim*Nb,dim))).T.reshape(-1)+(Nm-Nmi)
            aS = np.dot(np.dot(AP[i][:,col],D_array[i]),M_int[i]) # aS= (Id-D*T)\D*AP with M_Interaction=(Id-D*T)**-1. However, all have been transposed for sake of convenience
            # Save amplitude coefficients
            if self.cylamplitude:
                self.aS.append(aS)
                self.AP.append(AP[i][:,col])
            # Calculation of the excitation force for the entire array
            Fex[i,:,:] = np.dot(AP[i][:,col]+np.dot(aS,T[i]),G_array[i]) # Fex= G*aI with aI=AP+T*aS so the overall incident wave for the scattering problem
        self.Fex= Fex

    def Radiation(self, coord, Evanescent):
        """
        """
        Nb = len2(coord)
        Nm = self.Body.order.max()
        dof = self.Body.dof
        AR_iso = self.Body.AR
        Madd_iso = self.Body.Madd
        Crad_iso = self.Body.Crad
        k0 = self.wnumber
        freq = 2*pi/self.period
        #
        T = self.Transformation(coord, Nm, Evanescent)
        T, D_array, G_array, M_int = self.Interaction(T, Nm, 'R', Evanescent)
        # From now on we work for each wave-frequency
        Madd = np.zeros((len2(k0), Nb*dof, Nb*dof))
        Crad = np.zeros((len2(k0), Nb*dof, Nb*dof))
        # For later usage D, G, T and M_inter are transposed
        # Save amplitude coefficients
        if self.cylamplitude:
            self.aR= []
            self.AR= []
        for k in range(len2(k0)):
            Nmk = (len2(T[k])//Nb-1)//2
            dimk = 2*Nmk+1
            colT = np.array(np.arange(0,dimk*Nb,dimk).tolist()*dimk,dtype=int)+np.linspace(0,dimk,dimk*Nb,endpoint=False,dtype=int)
            colAR = (np.array([range(dimk)])+(Nm-Nmk)).reshape(-1)
            T_re = np.reshape(T[k][colT],(dimk,dimk*Nb**2)) # remember that T is already transposed so we are already getting T columns
            AR = np.dot(AR_iso[k][:,colAR],T_re) # each row of aR is constant dof.
            AR = np.reshape(AR,(dof,Nb,Nb*dimk))
            Fex_rad = np.zeros((dof,Nb,dof*Nb),dtype=complex)
            # Save amplitude coefficients
            if self.cylamplitude:
                    self.AR.append(AR)
                    aRaux = np.zeros(AR.shape,dtype=complex)
            for i in range(dof):
                aR = np.dot(np.dot(AR[i,:,:],D_array[k]),M_int[k])
                # Save amplitude coefficients
                if self.cylamplitude:
                    aRaux[i] = aR
                Fex_rad[i,:,:] = np.dot(AR[i,:,:]+np.dot(aR,T[k]),G_array[k])
            # Save amplitude coefficients
            if self.cylamplitude:
                self.aR.append(aRaux)
            Madd_t = Madd_iso[k,:,:].T.repeat(Nb,axis=0).reshape((dof,1,dof*Nb))
            Crad_t = Crad_iso[k,:,:].T.repeat(Nb,axis=0).reshape((dof,1,dof*Nb))
            Ones = np.eye(Nb,dtype=int).repeat(dof,axis=1) # will be used to distribute Crad_t over eye() type
            if self.Body.convention == 'N' :
                Frad_k = -(-freq[k]**2*Ones*Madd_t-1j*freq[k]*Ones*Crad_t)+Fex_rad
            else :
                Frad_k = -(-freq[k]**2*Ones*Madd_t+1j*freq[k]*Ones*Crad_t)+Fex_rad
            # Redistribution to a conventional way
            Frad_k = Frad_k.reshape((Nb*dof,Nb*dof), order= 'F').T
            Madd[k,:,:] = 1/freq[k]**2*Frad_k.real
            if self.Body.convention == 'N' :
                Crad[k,:,:] = 1/freq[k]*Frad_k.imag
            else :
                Crad[k,:,:] = -1/freq[k]*Frad_k.imag
        self.Madd= Madd
        self.Crad= Crad

    def Interaction(self, T, Nm, problem, Evanescent):
        """
        """
        dof = self.Body.dof
        dim = 2*Nm+1
        kl=self.wnumber_E
        k0 = self.wnumber
        if Evanescent : 
            Nb = T.shape[1]//(dim*(len(kl[0])+1))
        else :
            Nb = T.shape[1]//dim
        D = self.Body.D
        G = self.Body.G
        if problem == 'R' : # radiation
            TruncOrder = self.Body.truncorder[:,1]
        else : # scattering
            TruncOrder = self.Body.truncorder[:,0]
        # Interaction
        if Evanescent : 
            dim_e=dim*(len(kl[0])+1) # evanescent waves
        else : 
            dim_e=dim # progressive waves
        D_array = np.zeros((len2(k0),Nb*dim_e,Nb*dim_e),dtype=complex)
        G_array = np.zeros((len2(k0),Nb*dim_e,Nb*dof),dtype=complex) # G is already transposed, numb columns= Nb*dof ok!
        for i in range(Nb):
            D_array[:,dim_e*i:dim_e*(i+1),dim_e*i:dim_e*(i+1)] = D # diffraction transfer matrix for the isolated device
            G_array[:,dim_e*i:dim_e*(i+1),dof*i:dof*(i+1)] = G

        dimp = 2*TruncOrder+1
        T_trunc, D_trunc, G_trunc, M_inter = [], [], [], []

        for i in range(len2(k0)):
            rowcol = (np.array([range(dimp[i])]*Nb).transpose()+np.array(range(0,dim*Nb,dim))).transpose().reshape(-1)
            row, col = np.meshgrid(rowcol, rowcol, indexing='ij')
            ii = Nm-TruncOrder[i]
            D_trunc.append(D_array[i,row+ii,col+ii])
            G_trunc.append(G_array[i,rowcol+ii,:])
            T_trunc.append(T[i,row,col])
            if Evanescent : 
                M_inter.append(solve(np.eye(Nb*dim_e)-np.dot(T[i],D_array[i]),np.eye(Nb*dim_e)))
            else : 
                M_inter.append(solve(np.eye(Nb*dimp[i])-np.dot(T_trunc[i],D_trunc[i]),np.eye(Nb*dimp[i])))
        
        if Evanescent:
            return T, D_array, G_array, M_inter     # sans réduction
        else : 
            return T_trunc, D_trunc, G_trunc, M_inter # réduction 

    def Transformation(self, coord, Nm, Evanescent) :
        """
        """
        ##
        (k0, Nb, Nf, dim) = (self.wnumber, len2(coord), len2(self.wnumber), 2*Nm+1)
        if Evanescent : 
            kl=self.wnumber_E
        # Geometric properties
        (row, col) = np.meshgrid(range(Nb), range(Nb), indexing = 'ij')
        Dx = (coord[:,0]*np.ones((Nb,Nb))).T-coord[:,0]*np.ones((Nb,Nb))
        Dy = (coord[:,1]*np.ones((Nb,Nb))).T-coord[:,1]*np.ones((Nb,Nb))
        L = np.sqrt(Dx**2+Dy**2)[row<col]
        alf = np.arctan2(Dy,Dx)[row<col]
        ij_i = row[row<col]
        ij_j = col[row<col]
        ##
        (D, kk0, nu) = np.meshgrid(L, k0, range(dim), indexing='ij')
        H = jv(nu, kk0*D)-1j*yv(nu, kk0*D)
        if self.Body.convention == 'N' :
            H = np.conj(H)
        del(D, kk0, nu)
        if Evanescent : 
            (D, kkl, nu) = np.meshgrid(L, kl, range(dim), indexing='ij')
            Ke = kv(nu, kkl*D)
            del(D, kkl, nu)
        
        # Transformation matrix
        (p, q) = np.meshgrid(range(-Nm, Nm+1), range(-Nm, Nm+1), indexing = 'ij')
        mapk0 = np.repeat(range(Nf), dim*dim).reshape((Nf, dim, dim))
        sign = (np.tri(dim).T-np.tri(dim))**(q-p)
        if Evanescent : 
            dim_e=dim*(len(kl[0])+1) # evanescent waves
        else : 
            dim_e=dim # progressive waves
        T = np.zeros((Nf, Nb*dim_e, Nb*dim_e),dtype=complex)
        for i in range(len2(L)):
            r = dim*ij_i[i]
            c = dim*ij_j[i]
            # Assembling Tij and Tji
            Tij= sign*H[i][mapk0, abs(q-p)]*np.exp((q-p)*alf[i]*1j)
            T[:, r:r+dim, c:c+dim]= Tij # Tij
            T[:, c:c+dim, r:r+dim]= Tij*np.exp((q-p)*pi*1j) # Tji

            if Evanescent : 
                for j in range(len2(kl[0])) :
                    mode_e=(j+1)*Nb*dim
                    Tij_e= sign*Ke[i, j, abs(q-p)]*np.exp((q-p)*alf[i]*1j)*(-1.0)**p
                    T[:, mode_e+r:mode_e+r+dim, mode_e+c:mode_e+c+dim]= Tij_e # Tij
                    T[:, mode_e+c:mode_e+c+dim, mode_e+r:mode_e+r+dim]= Tij_e*np.exp((q-p)*pi*1j) # Tji

        #Transpose for later usage since we finally  decided
        #to work with transposed(T). D and G for the isolated device
        #have been already obtained transposed as well as
        #AP, aS, AR and aR
        T = np.transpose(T,(0,2,1))
        return T
    
    def RAO(self, directory) : 
        Madd=self.Madd
        Crad=self.Crad
        Fex=self.Fex
        freq = 2*pi/self.period
        Nfreq = freq.shape[0]
        Ndir = Fex.shape[1]  
        Ndof = Fex.shape[2]  # 6*n bodies 
        Mass=0
        if len(glob(os.path.join(directory,'Inertia.dat'))) > 0 :
            with open(os.path.join(directory,'Inertia.dat'), 'r') as f_Mass :
                Mass = np.loadtxt(f_Mass)
        Kh=0
        if len(glob(os.path.join(directory,'Kh.dat'))) > 0 :
            with open(os.path.join(directory,'Kh.dat'), 'r') as f_Kh :
                Kh = np.loadtxt(f_Kh)
        RAO = np.zeros((Nfreq, Ndir, Ndof), dtype=complex)  

        for i in range(Nfreq):
            mat_A = -(Mass + Madd[i,:,:]) * freq[i]**2 - 1j * freq[i] * Crad[i,:,:] + Kh  
            for beta in range(Ndir) : 
                RAO[i,beta, :] = np.linalg.lstsq(mat_A, Fex[i, beta, :], rcond=None)[0]

        self.RAO= RAO

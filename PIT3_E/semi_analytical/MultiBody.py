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
                                    # AP[w, beta, body, (2*Nm+1)*(l+1)+m] = np.exp(1j*mode[w, beta, body, m]*(pi/2-matdir[w, beta, body, m]))
                                    # AP[w, beta, body, (2*Nm+1)*(l+1)+m] *= np.exp(1j*kl[w][l]*(X[w, beta, body, m]*np.cos(matdir[w, beta, body, m])+Y[w, beta, body, m]*np.sin(matdir[w, beta, body, m])))
                                    AP[w, beta, body, (2*Nm+1)*(l+1)+m] = 0
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
            if Evanescent : 
                Ne=len2(kl[0])
            else :
                Ne=0
            Nmi = (len2(T[i])//(Nb*(Ne+1))-1)//2
            dimi = 2*Nmi+1
            col = (np.array([range(dimi*(Ne+1))]*Nb).T+np.array(range(0,dim*Nb*(Ne+1),dim*(Ne+1)))).T.reshape(-1)+(Nm-Nmi)*(Ne+1)
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
        kl = self.wnumber_E
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
            if Evanescent : 
                Ne=len2(kl[0])
            else :
                Ne=0
            Nmk = (len2(T[k])//(Nb*(Ne+1))-1)//2
            dimk = 2*Nmk+1
            colT = np.array(np.arange(0,dimk*(Ne+1)*Nb,dimk*(Ne+1)).tolist()*dimk*(Ne+1),dtype=int)+np.linspace(0,dimk*(Ne+1),dimk*(Ne+1)*Nb,endpoint=False,dtype=int)
            colAR = (np.array([range(dimk*(Ne+1))])+(Nm-Nmk)*(Ne+1)).reshape(-1)
            # print("T[k].shape:", T[k].shape)
            # print("colT.shape:", colT.shape)
            # print("colT[:10]:", colT[:10])
            # print("colAR[:10]:", colAR[:10])
            T_re = np.reshape(T[k][colT],(dimk*(Ne+1),dimk*(Ne+1)*Nb**2)) # remember that T is already transposed so we are already getting T columns
            AR = np.dot(AR_iso[k][:,colAR],T_re) # each row of aR is constant dof.
            AR = np.reshape(AR,(dof,Nb,Nb*dimk*(Ne+1)))
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
            # print("dim_e*i:dim_e*(i+1),dim_e*i:dim_e*(i+1)", dim_e*i,dim_e*(i+1),dim_e*i,dim_e*(i+1))

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
        # print(H.shape)

        if Evanescent : 
            kl=self.wnumber_E
            Ne=len2(kl[0])
            dim_e=dim*(Ne+1) # evanescent waves
        
            # METHOD 1 
            Ke=np.zeros((Ne, len2(L), Nf, dim)) 
            for E in range(Ne) : 
                (D, kkl, nu) = np.meshgrid(L, kl[:,E], range(dim), indexing='ij') 
                Ke[E]=kv(nu, kkl*D)
                del(D, kkl, nu)
                # print(Ke.shape)

             # METHOD 2
            # D = L[:, None, None]             # shape (Nd, 1, 1)
            # kkl = kl[None, :, :]             # shape (1, Nf, Ne)
            # nu = np.arange(dim)[None, None, :, None]  # shape (1, 1, dim, 1)
            # Dkkl = (D * kkl)[:, :, None, :]            # shape (Nd, Nf, 1, Ne)
            # Ke_full = kv(nu, Dkkl)  # shape (Nd, Nf, dim, Ne)
            # # Ke = Ke_full.reshape(len2(L), Nf, dim * Ne) 

            # pe = np.tile(p, (Ne, Ne))  # Répète p Ne fois verticalement et horizontalement
            # qe = np.tile(q, (Ne, Ne))

            #  METHOD 2
            # idx = np.arange(dim * Ne)
            # row_idx = np.repeat(idx, dim * Ne)
            # col_idx = np.tile(idx, dim * Ne)
            # pe = row_idx % dim
            # qe = col_idx % dim
            # le1 = row_idx // dim
            # le2 = col_idx // dim
            # mask = (le1 == le2)
            # # sign_e = (np.tri(dim*Ne).T-np.tri(dim*Ne))**(qe-pe)
            # sign_e = np.where(qe >= pe, 1.0, -1.0)
        else : 
            dim_e=dim # progressive waves

        # Transformation matrix
        (p, q) = np.meshgrid(range(-Nm, Nm+1), range(-Nm, Nm+1), indexing = 'ij')
        mapk0 = np.repeat(range(Nf), dim*dim).reshape((Nf, dim, dim))
        sign = (np.tri(dim).T-np.tri(dim))**(q-p)            
            
        T = np.zeros((Nf, Nb*dim_e, Nb*dim_e),dtype=complex)
        for i in range(len2(L)):
            r = dim_e*ij_i[i]
            c = dim_e*ij_j[i]
            Tij= sign*H[i][mapk0, abs(q-p)]*np.exp((q-p)*alf[i]*1j)
            # print(Tij.shape)
            T[:, r:r+dim, c:c+dim]= Tij # Tij
            T[:, c:c+dim, r:r+dim]= Tij*np.exp((q-p)*pi*1j) # Tji
            # print("r:r+dim, c:c+dim", r,r + dim, c, c + dim)
            if Evanescent : 
                # METHOD 1
                for E in range(Ne) :
                    r_e=dim_e*ij_i[i]+(E+1)*dim
                    c_e=dim_e*ij_j[i]+(E+1)*dim
                    Tij_e = sign * Ke[E][i][mapk0,  abs(q-p)]  * np.exp((q - p)*alf[i]*1j) * (-1.0)**p
                    T[:, r_e:r_e + dim, c_e:c_e + dim] = Tij_e
                    T[:, c_e:c_e + dim, r_e:r_e + dim] = Tij_e * np.exp((q - p)*np.pi*1j)  
                    # print("Tij_e", Tij_e.shape)
                    # print("ij", ij_i[i], ij_j[i])
                    # print("E=", E)
                    # print("r_e:r_e + dim, c_e:c_e + dim", r_e,r_e + dim, c_e, c_e + dim)
                
                    # print("Tij_e", Tij_e.shape, Tij_e[-1, :, :])
                    # print("formule", kv(0, kl[-1,E]*L[i])*(-1))
                    # print("formule", kv(1, kl[-1,E]*L[i])*(-1) * np.exp(alf[i]*1j))
                    # print("formule", kv(4, kl[-1,E]*L[i])*(-1) * np.exp(4*alf[i]*1j))

                # METHOD 2
                # r_e = dim_e * ij_i[i] + dim
                # c_e = dim_e * ij_j[i] + dim
                # for w in range(Nf):
                    # Ke_w = Ke_full[i, w]  # shape (dim, Ne)
                    # # Tij_e = sign_e * Ke_w[np.abs(qe-pe), le1]  * np.exp((qe - pe)*alf[i]*1j) * (-1.0)**pe
                    # Tij_e_flat = sign_e * Ke_w[np.abs(qe - pe), le1] \
                    #             * np.exp((qe - pe) * alf[i] * 1j) \
                    #             * (-1.0) ** pe
                    # Tij_e = np.zeros((dim * Ne, dim * Ne), dtype=complex)
                    # Tij_e[row_idx, col_idx] = Tij_e_flat.flatten()
                    # # Tij_e[row_idx[mask], col_idx[mask]] = (sign_e 
                    # # * Ke_w[np.abs(qe[mask] - pe[mask]), le1[mask]] 
                    # # * np.exp((qe[mask] - pe[mask]) * alf[i] * 1j) * (-1.0) ** pe[mask])
                    # T[w, r_e:r_e + dim*Ne, c_e:c_e + dim*Ne] = Tij_e
                    # T[w, c_e:c_e + dim*Ne, r_e:r_e + dim*Ne] = Tij_e * np.exp((qe - pe).reshape(dim*Ne, dim*Ne)*np.pi*1j)  

                
                
        #Transpose for later usage since we finally  decided
        #to work with transposed(T). D and G for the isolated device
        #have been already obtained transposed as well as
        #AP, aS, AR and aR
        T = np.transpose(T,(0,2,1))
        # print(T.shape)
        # print(T[0])
        # Écriture dans le fichier
        # with open("T.dat", 'w') as f:
        #     f.write("# Matrice complexe (ligne par ligne)\n")
        #     for row in T[0]:
        #         line = "\t".join(f"{val.real:.6f}+{val.imag:.6f}j" for val in row)
        #         f.write(line + "\n")
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

    def Kochin(self, theta) : 
        aR=self.aR
        aS=self.aS
        RAO=self.RAO
        # print("Radiation")
        # print(len(self.aR), len(self.aR[0]), len(self.aR[0][0]), len(self.aR[0][0][0]), len(self.aR[1][0][0]))
        # print("Scattering")
        # print(len(self.aS), len(self.aS[0]), len(self.aS[0][0]), len(self.aS[1][0]))
        
        freq = 2*pi/self.period
        Nfreq = freq.shape[0]
        Ndir = len(self.aS[0])
        Ndof = len(self.aR[0])
        Nb=len(self.aR[0][0]) # 6*n bodies 
        KochinR = np.zeros((Nfreq, Ndof*Nb, len(theta)), dtype=complex)
        KochinS = np.zeros((Nfreq, Ndir, len(theta)), dtype=complex)
        g=9.81
        for ind, w in enumerate(freq):
            NmR = len(aR[ind][0][0])
            NmS = len(aS[ind][0])
            M_R=int((NmR-1)/2)
            M_S=int((NmS-1)/2)
            for k, angle in enumerate(theta) :
                for i, mode in enumerate(range(-M_R, M_R+1)):
                    for body in range(Nb) :
                        for dof in range(Ndof) :
                            index_dof = body * Ndof + dof
                            KochinR[ind, index_dof, k]+= (-1)**mode *(1j*w/g) *aR[ind][dof][body][i]*np.exp(1j*mode*angle)
                            # KochinR[ind, index_dof, k]+= 1j*w/g *aR[ind][dof][body][i]*np.exp(1j*mode*angle)*RAO[ind,0,body]*np.pi/180
                for i, mode in enumerate(range(-M_S, M_S+1)):
                    for beta in range(Ndir) :
                        KochinS[ind, beta, k]+= (-1)**mode *1j*w/g *aS[ind][beta][i]*np.exp(1j*mode*angle)
        self.KochinR=KochinR
        self.KochinS=KochinS


    def FreeSurface(self, coord, Nx, Ny, Lx, Ly) : 
        aR=self.aR
        aS=self.aS
        k0=self.wnumber
        direction = self.dir
        freq = 2 * pi / self.period
        Nfreq = freq.shape[0]
        Ndir = len(self.aS[0])
        Ndof = len(self.aR[0])
        Nb = len(self.aR[0][0]) 
        ETA = np.zeros((Nfreq, Ndir, Nx, Ny), dtype=complex)
        ETA_R = np.zeros((Nfreq, Ndof*Nb, Nx, Ny), dtype=complex)
        phiR = np.zeros((Nfreq, Nx, Ny), dtype=complex)
        phiS = np.zeros((Nfreq, Ndir, Nx, Ny), dtype=complex)
        ETA_S = np.zeros_like(phiS)
        g = 9.81

        x = np.linspace(-Lx/2, Lx/2, Nx)
        y = np.linspace(-Ly/2, Ly/2, Ny)
        X, Y = np.meshgrid(x, y, indexing='ij')  # (Nx, Ny)

        for ind, w in enumerate(freq):
            NmR = len(aR[ind][0][0])
            NmS = len(aS[ind][0])
            M_R = int((NmR - 1) / 2)
            M_S = int((NmS - 1) / 2)
            k = k0[ind]
            for body in range(Nb):
                dx = X - coord[body, 0]
                dy = Y - coord[body, 1]
                L = np.sqrt(dx**2 + dy**2)
                alpha = np.arctan2(dy, dx)
                for dof in range(Ndof):
                    index = body * Ndof + dof
                    for idx_mode, mode in enumerate(range(-M_R, M_R + 1)):
                        Hm = jv(mode, k * L) + 1j * yv(mode, k * L)  # Hankel function type 1
                        contribution = aR[ind][dof][body][idx_mode] * Hm * np.exp(1j * mode * alpha)
                        ETA_R[ind, index] += contribution
                        phiR[ind] += contribution
                for beta in range(Ndir):
                    for idx_mode, mode in enumerate(range(-M_S, M_S + 1)):
                        Hm = jv(mode, k * L) + 1j * yv(mode, k * L)
                        phiS[ind, beta] += aS[ind][beta][idx_mode] * Hm * np.exp(1j * mode * alpha)

            for beta in range(Ndir):
                phiI = np.exp(1j * k * (X * np.cos(direction[beta]) + Y * np.sin(direction[beta])))
                ETA[ind, beta] = (phiS[ind, beta] + phiR[ind]+phiI) * (1j * w / g)
            ETA_R[ind] *= (1j * w / g)
            ETA_S[ind] = phiS[ind]* (1j * w / g)
        self.ETA=ETA
        self.ETA_R=ETA_R 
        self.ETA_S=ETA_S
        
        # for ind, w in enumerate(freq):
        #     phiR = np.zeros((Nx, Ny), dtype=complex)
        #     phiS = np.zeros((Ndir, Nx, Ny), dtype=complex)
        #     NmR = len(aR[ind][0][0])
        #     NmS = len(aS[ind][0])
        #     M_R=int((NmR-1)/2)
        #     M_S=int((NmS-1)/2)
        #     for i in range(Nx):
        #         for j in range(Ny):
        #             for body in range(Nb) :
        #                 L=np.sqrt((x[i]-coord[body, 1])**2+(y[j]-coord[body, 2])**2)
        #                 alpha=np.arctan2((y[j]-coord[body, 2]), (x[i]-coord[body, 1]))
        #                 for dof in range(Ndof) :
        #                     for k, mode in enumerate(range(-M_R, M_R+1)):
        #                         phiR[i, j] = phiR[i, j] + aR[ind][dof][body][k]*(jv(mode, k0[ind]*L)+yv(mode,k0[ind]*L))*np.exp(1j*mode*alpha)
        #                 for beta in range(Ndir) :
        #                     for k, mode in enumerate(range(-M_S, M_S+1)):
        #                         phiS[beta, i, j] = phiS[beta, i, j] + aS[ind][beta][k]*(jv(mode, k0[ind]*L)+yv(mode,k0[ind]*L))*np.exp(1j*mode*alpha)
        #             for beta in range(Ndir) :
        #                 phiI=np.exp(1j*k0[ind]*(x[i]*np.cos(direction[beta])+y[j]*np.sin(direction[beta])))
        #                 ETA[ind, beta, i,j] = phiS[beta, i, j]+phiR[i, j]+phiI
        #                 ETA[ind, beta, i,j] *= 1j*w/g
                    

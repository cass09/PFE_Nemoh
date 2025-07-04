"""
"""

# Imports
from math import pi
import numpy as np
from numpy.linalg import solve
from scipy.special import jv, yv, kv
from glob import glob
import os
import cmath
from capytaine.bem.cylindrical_waves import WNumber_E
# Settings


class MultiBody(object):

    def __init__(self, directions, iBody, Ndof, Madd_iso, Crad_iso, cylamplitude= True):
        """
        """
        # cylamplitude = True ==> will save cylindrical amplitude coefficients
        self.cylamplitude = cylamplitude
        # wave data from iBody
        self.depth = iBody.depth
        self.period = iBody.period
        self.wnumber = iBody.wnumber
        self.L=iBody.Nmodes_E
        self.wnumber_E = WNumber_E(self.period, self.depth, self.L, One=False)
        # wave directions
        self.dir = directions
        # Body instance
        self.Body = iBody
        self.dof = Ndof
        self.Madd_iso = Madd_iso
        self.Crad_iso = Crad_iso

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

    def Scattering(self, coord, L):
        """
        """
        Nb = len(coord)
        Nm = self.Body.order.max()
        dim = 2*Nm+1
        k0 = self.wnumber
        # kl = self.wnumber_E
        direction = self.dir
        #
        T = self.Transformation(coord, Nm, L)
        T, D_array, G_array, M_int = self.Interaction(T, Nm, 'S', L)
        # Get the ambient planar wave, AP, for the scattering problem for the entire array
        matk0,matdir,X,mode = np.meshgrid(k0, direction, coord[:,0], range(-Nm,Nm+1),indexing='ij')
        Y = np.meshgrid(k0, direction, coord[:,1], range(-Nm,Nm+1), indexing='ij')[2]
        # Calculate AP using (k0,dir,bodycoord,mode)
        if L>0 :
            Ntot=(2*Nm+1)*(L+1)
            AP = np.zeros((len(k0),len(direction),len(coord), Ntot), dtype=complex)
            for w in range(len(k0)) : 
                for beta in range(len(direction)):
                    for body in range(len(coord)):
                        for m in range(dim):
                            AP[w, beta, body, m] = np.exp(1j*mode[w, beta, body, m]*(pi/2-matdir[w, beta, body, m]))
                            AP[w, beta, body, m] *= np.exp(1j*matk0[w, beta, body, m]*(X[w, beta, body, m]*np.cos(matdir[w, beta, body, m])+Y[w, beta, body, m]*np.sin(matdir[w, beta, body, m])))
                            for l in range(L):
                                # AP[w, beta, body, (2*Nm+1)*(l+1)+m] = np.exp(1j*mode[w, beta, body, m]*(pi/2-matdir[w, beta, body, m]))
                                # AP[w, beta, body, (2*Nm+1)*(l+1)+m] *= np.exp(1j*kl[w][l]*(X[w, beta, body, m]*np.cos(matdir[w, beta, body, m])+Y[w, beta, body, m]*np.sin(matdir[w, beta, body, m])))
                                AP[w, beta, body, (2*Nm+1)*(l+1)+m] = 0
        else : 
            AP = np.exp(1j*mode*(pi/2-matdir))
            AP *= np.exp(1j*matk0*(X*np.cos(matdir)+Y*np.sin(matdir)))
        # else :
        #     AP = np.exp(-1j*mode*(pi/2+matdir))
        #     AP *= np.exp(-1j*matk0*(X*np.cos(matdir)+Y*np.sin(matdir)))
        AP = np.reshape(AP,(len(k0),len(direction),(2*Nm+1)*Nb*(L+1)))
        Fex = np.zeros((len(k0),len(direction),len(G_array[0][0])), dtype=complex)
        # Save amplitude coefficients
        if self.cylamplitude:
            self.aS = []
            self.AP = []
        for i in range(len(k0)):
            Nmi = (len(T[i])//(Nb*(L+1))-1)//2
            dimi = 2*Nmi+1
            col = (np.array([range(dimi*(L+1))]*Nb).T+np.array(range(0,dim*Nb*(L+1),dim*(L+1)))).T.reshape(-1)+(Nm-Nmi)*(L+1)
            aS = np.dot(np.dot(AP[i][:,col],D_array[i]),M_int[i]) # aS= (Id-D*T)\D*AP with M_Interaction=(Id-D*T)**-1. However, all have been transposed for sake of convenience
            # Save amplitude coefficients
            if self.cylamplitude:
                self.aS.append(aS)
                self.AP.append(AP[i][:,col])
            # Calculation of the excitation force for the entire array
            Fex[i,:,:] = np.dot(AP[i][:,col]+np.dot(aS,T[i]),G_array[i]) # Fex= G*aI with aI=AP+T*aS so the overall incident wave for the scattering problem
        self.Fex= Fex
        
        
    def Radiation(self, coord, L):
        """
        """
        Nb = len(coord)
        Nm = self.Body.order.max()
        dof = self.dof
        AR_iso = self.Body.AR
        Madd_iso = self.Madd_iso
        Crad_iso = self.Crad_iso
        k0 = self.wnumber
    
        # kl = self.wnumber_E
        freq = 2*pi/self.period
        #
        T = self.Transformation(coord, Nm, L)
        T, D_array, G_array, M_int = self.Interaction(T, Nm, 'R', L)
        # From now on we work for each wave-frequency
        Madd = np.zeros((len(k0), Nb*dof, Nb*dof))
        Crad = np.zeros((len(k0), Nb*dof, Nb*dof))
        # For later usage D, G, T and M_inter are transposed
        # Save amplitude coefficients
        if self.cylamplitude:
            self.aR= []
            self.AR= []
        for k in range(len(k0)):
            Nmk = (len(T[k])//(Nb*(L+1))-1)//2
            dimk = 2*Nmk+1
            colT = np.array(np.arange(0,dimk*(L+1)*Nb,dimk*(L+1)).tolist()*dimk*(L+1),dtype=int)+np.linspace(0,dimk*(L+1),dimk*(L+1)*Nb,endpoint=False,dtype=int)
            colAR = (np.array([range(dimk*(L+1))])+(Nm-Nmk)*(L+1)).reshape(-1)
            # print("T[k].shape:", T[k].shape)
            # print("colT.shape:", colT.shape)
            # print("colT[:10]:", colT[:10])
            # print("colAR[:10]:", colAR[:10])
            T_re = np.reshape(T[k][colT],(dimk*(L+1),dimk*(L+1)*Nb**2)) # remember that T is already transposed so we are already getting T columns
            # print("T_re.shape:", T_re.shape)
            # print("T", T[k])
            # print("T_re", T_re)
            AR = np.dot(AR_iso[k][:,colAR],T_re) # each row of aR is constant dof.
            # print("shape", AR.shape, AR_iso[k].shape)
            AR = np.reshape(AR,(dof,Nb,Nb*dimk*(L+1)))
            # print("AR.shape", AR.shape)
            Fex_rad = np.zeros((dof,Nb,dof*Nb),dtype=complex)
            # Save amplitude coefficients
            if self.cylamplitude:
                    self.AR.append(AR)
                    aRaux = np.zeros(AR.shape,dtype=complex)
            for i in range(dof):
                aR = np.dot(np.dot(AR[i,:,:],D_array[k]),M_int[k])
                # print("aR.shape", aR.shape)
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
            # if self.Body.convention == 'N' :
            Frad_k = -(-freq[k]**2*Ones*Madd_t-1j*freq[k]*Ones*Crad_t)+Fex_rad
            # else :
            #     Frad_k = -(-freq[k]**2*Ones*Madd_t+1j*freq[k]*Ones*Crad_t)+Fex_rad
            # Redistribution to a conventional way
            Frad_k = Frad_k.reshape((Nb*dof,Nb*dof), order= 'F').T
            Madd[k,:,:] = 1/freq[k]**2*Frad_k.real
            # if self.Body.convention == 'N' :
            Crad[k,:,:] = 1/freq[k]*Frad_k.imag
            # else :
            #     Crad[k,:,:] = -1/freq[k]*Frad_k.imag
        self.Madd= Madd
        self.Crad= Crad

    def Interaction(self, T, Nm, problem, L):
        """
        """
        dof = self.dof
        dim = 2*Nm+1
        # kl=self.wnumber_E
        k0 = self.wnumber
        
        Nb = T.shape[1]//(dim*(L+1))
       
        D = self.Body.D
        G = self.Body.G
        if problem == 'R' : # radiation
            TruncOrder = self.Body.truncorder[:,1]
        else : # scattering
            TruncOrder = self.Body.truncorder[:,0]
        # Interaction
        dim_e=dim*(L+1) # evanescent waves
      
        D_array = np.zeros((len(k0),Nb*dim_e,Nb*dim_e),dtype=complex)
        G_array = np.zeros((len(k0),Nb*dim_e,Nb*dof),dtype=complex) # G is already transposed, numb columns= Nb*dof ok!
        for i in range(Nb):
            D_array[:,dim_e*i:dim_e*(i+1),dim_e*i:dim_e*(i+1)] = D # diffraction transfer matrix for the isolated device
            G_array[:,dim_e*i:dim_e*(i+1),dof*i:dof*(i+1)] = G
            # print("dim_e*i:dim_e*(i+1),dim_e*i:dim_e*(i+1)", dim_e*i,dim_e*(i+1),dim_e*i,dim_e*(i+1))

        dimp = 2*TruncOrder+1
        T_trunc, D_trunc, G_trunc, M_inter = [], [], [], []

        for i in range(len(k0)):
            rowcol = (np.array([range(dimp[i])]*Nb).transpose()+np.array(range(0,dim*Nb,dim))).transpose().reshape(-1)
            row, col = np.meshgrid(rowcol, rowcol, indexing='ij')
            ii = Nm-TruncOrder[i]
            D_trunc.append(D_array[i,row+ii,col+ii])
            G_trunc.append(G_array[i,rowcol+ii,:])
            T_trunc.append(T[i,row,col])
            if L>0 : 
                M_inter.append(solve(np.eye(Nb*dim_e)-np.dot(T[i],D_array[i]),np.eye(Nb*dim_e)))
            else : 
                M_inter.append(solve(np.eye(Nb*dimp[i])-np.dot(T_trunc[i],D_trunc[i]),np.eye(Nb*dimp[i])))
        
        if L>0:
            return T, D_array, G_array, M_inter     # sans réduction
        else : 
            return T_trunc, D_trunc, G_trunc, M_inter # réduction 

    def Transformation(self, coord, Nm, L) :
        """
        """
        ##
        (k0, Nb, Nf, dim) = (self.wnumber, len(coord), len(self.wnumber), 2*Nm+1)  
        # Geometric properties
        (row, col) = np.meshgrid(range(Nb), range(Nb), indexing = 'ij')
        Dx = (coord[:,0]*np.ones((Nb,Nb))).T-coord[:,0]*np.ones((Nb,Nb))
        Dy = (coord[:,1]*np.ones((Nb,Nb))).T-coord[:,1]*np.ones((Nb,Nb))
        Lenght = np.sqrt(Dx**2+Dy**2)[row<col]
        alf = np.arctan2(Dy,Dx)[row<col]
        ij_i = row[row<col]
        ij_j = col[row<col]
        ##
        (D, kk0, nu) = np.meshgrid(Lenght, k0, range(dim), indexing='ij')
        H = jv(nu, kk0*D)-1j*yv(nu, kk0*D)
        # if self.Body.convention == 'N' :
        H = np.conj(H)
        del(D, kk0, nu)
        # print(H.shape)

        if L>0 : 
            kl=self.wnumber_E
            dim_e=dim*(L+1) # evanescent waves
        
            # METHOD 1 
            Ke=np.zeros((L, len(Lenght), Nf, dim)) 
            for E in range(L) : 
                (D, kkl, nu) = np.meshgrid(Lenght, kl[:,E], range(dim), indexing='ij') 
                Ke[E]=kv(nu, kkl*D)
                del(D, kkl, nu)
                # print(Ke.shape)

             # METHOD 2
            # D = L[:, None, None]             # shape (Nd, 1, 1)
            # kkl = kl[None, :, :]             # shape (1, Nf, Ne)
            # nu = np.arange(dim)[None, None, :, None]  # shape (1, 1, dim, 1)
            # Dkkl = (D * kkl)[:, :, None, :]            # shape (Nd, Nf, 1, Ne)
            # Ke_full = kv(nu, Dkkl)  # shape (Nd, Nf, dim, Ne)
            # # Ke = Ke_full.reshape(len(L), Nf, dim * Ne) 

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
        for i in range(len(Lenght)):
            r = dim_e*ij_i[i]
            c = dim_e*ij_j[i]
            Tij= sign*H[i][mapk0, abs(q-p)]*np.exp((q-p)*alf[i]*1j)
            # print(Tij.shape)
            T[:, r:r+dim, c:c+dim]= Tij # Tij
            T[:, c:c+dim, r:r+dim]= Tij*np.exp((q-p)*pi*1j) # Tji
            # print("r:r+dim, c:c+dim", r,r + dim, c, c + dim)
            if L>0 : 
                # METHOD 1
                for E in range(L) :
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
        B_PTO=0
        K_PTO=0
        for i in range(Nfreq):
            mat_A = -(Mass + Madd[i,:,:]) * freq[i]**2 - 1j * freq[i] * (Crad[i,:,:]+B_PTO) + (Kh +K_PTO) 
            for beta in range(Ndir) : 
                RAO[i,beta, :] = np.linalg.lstsq(mat_A, Fex[i, beta, :], rcond=None)[0]

        self.RAO= RAO

    def Kochin(self, theta) : 
        aR=self.aR
        AR_iso = self.Body.AR
        aS=self.aS
        k0=self.wnumber
        d=self.depth
        # RAO=self.RAO
        # print("Radiation")
        # print(len(self.aR), len(self.aR[0]), len(self.aR[0][0]), len(self.aR[0][0][0]), len(self.aR[1][0][0]))
        # print("Scattering")
        # print(len(self.aS), len(self.aS[0]), len(self.aS[0][0]), len(self.aS[1][0]))
        freq = 2*pi/self.period
        Nfreq = freq.shape[0]
        Ndir = len(self.aS[0])
        Ndof = len(self.aR[0])
        Nb=len(self.aR[0][0])  
        KochinR = np.zeros((Nfreq, Ndof*Nb, len(theta)), dtype=complex)
        KochinS = np.zeros((Nfreq, Ndir, len(theta)), dtype=complex)
        g=9.81
        for ind, w in enumerate(freq):
            NmR = len(aR[ind][0][0])
            NmS = len(aS[ind][0])
            M_R=int((NmR-1)/(2*Nb))
            M_S=int((NmS-1)/(2*Nb))
            k_factor=(0.5*np.exp(k0[ind]*d)/(np.sinh(k0[ind]*d)*k0[ind]))
            for k, angle in enumerate(theta) :
                for body in range(Nb) :
                    for dof in range(Ndof) :
                        index_dof = body * Ndof + dof
                        coef_iso=AR_iso[ind][dof][0:2*M_R+1] #/(1j*w)
                        for j in range(Nb):  # loop to sum the contribution
                            coef=aR[ind][dof][body][(2*M_R+1)*j:(2*M_R+1)*(j+1)]
                            for mode in range(-M_R, M_R + 1):
                                # contribution = coef[M_R+mode] * (-1j)**(mode)  * np.exp(1j * mode * angle)
                                if j==body : 
                                    contribution = coef_iso[M_R+mode] * (-1j)**(mode)  * np.exp(1j * mode * angle)
                                
                                # if mode>0:
                                #     contribution+=coef[M_R-mode] * (-1j)**(mode) * np.exp(-1j * mode * angle)*(-1)**(mode)
                                #     if j==body :
                                #         contribution+=coef_iso[M_R-mode] * (-1j)**(mode) * np.exp(-1j * mode * angle)*(-1)**(mode)
                                KochinR[ind, index_dof, k]+= contribution
                        KochinR[ind, index_dof, k]= KochinR[ind, index_dof, k]*(1j*w/g)*(-1j*g/w)/(1j*w) *k_factor
                            # KochinR[ind, index_dof, k]+= (-1j)**mode *(1j*w/g) *aR[ind][dof][body][i]*np.exp(1j*mode*angle)
                                # KochinR[ind, index_dof, k]+= (1j)**mode  *aR[ind][dof][body][idx_mode]*np.exp(1j*mode*angle)
                                # KochinR[ind, index_dof, k]+= (1j)**mode *(1j*w/g)  *np.conj(aR[ind][dof][body][idx_mode])*np.exp(-1j*mode*angle)
                for beta in range(Ndir) :
                    for j in range(Nb):  
                        coef=aS[ind][beta][(2*M_S+1)*j:(2*M_S+1)*(j+1)]
                        for mode in range(0, M_S + 1):
                            contribution = coef[M_S+mode] * (-1j)**(mode) * np.exp(1j * mode * angle)
                            if mode>0:
                                contribution+=coef[M_S-mode] * (-1j)**(-mode) * np.exp(-1j * mode * angle)*(-1)**(mode)
                            KochinS[ind, beta, k] += contribution 
                    KochinS[ind, beta, k]=KochinS[ind, beta, k]*(1j*w/g)*(1j*g/w) *k_factor
                            # idx_mode = j * (2*M_S+1) + (mode + M_S)  
                            # KochinS[ind, beta, k]+= (1j)**mode*(1j*w/g)  *np.conj(aS[ind][beta][idx_mode])*np.exp(-1j*mode*angle)
        scale=1
        # scale=np.sqrt(2/np.pi)*np.exp(1j*np.pi/4)
        # self.KochinR=np.sqrt(2/np.pi)*np.exp(1j*np.pi/4)*KochinR
        self.KochinR=KochinR*scale
        self.KochinS=KochinS*scale
        

    def FreeSurface(self, coord, Nx, Ny, Lx, Ly, Ne) : 
        aR=self.aR
        AR_iso = self.Body.AR
        # print("aR", len(aR), len(aR[0]), len(aR[0][0]), len(aR[0][0][0]))
        # print("aR iso", AR_iso.shape)
        aS=self.aS
        k0=self.wnumber
        kl=self.wnumber_E
        depth=self.depth
        direction = self.dir
        freq = 2 * pi / self.period
        Nfreq = freq.shape[0]
        Ndir = len(aS[0])
        Ndof = len(aR[0])
        Nb = len(aR[0][0]) 
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
            print("w=", w)
            # number of progressive modes M (max=5)
            # can be lower for some frequencies
            NmR = len(aR[ind][0][0])
            NmS = len(aS[ind][0])
            M_R = int((NmR/(Nb*(Ne+1)) - 1) / 2)
            M_S = int((NmS/(Nb*(Ne+1)) - 1) / 2)
            k = k0[ind]
            kll=kl[ind,:]
            # loop on the bodies
            for body in range(Nb):
                # for radiation problem
                dx_i = X - coord[body, 0]
                dy_i = Y - coord[body, 1]
                L_i = np.sqrt(dx_i**2 + dy_i**2)
                alpha_i = np.arctan2(dy_i, dx_i)
                for dof in range(Ndof):
                    index = body * Ndof + dof
                    coef_iso=AR_iso[ind][dof][0:2*M_R+1] #/(1j*w)
                    
                    for j in range(Nb):  # loop to sum the contribution
                        dx = X - coord[j, 0]
                        dy = Y - coord[j, 1]
                        L = np.sqrt(dx**2 + dy**2)
                        alpha = np.arctan2(dy, dx)
                        coef=aR[ind][dof][body][(2*M_R+1)*j*(Ne+1):(2*M_R+1)*(j*(Ne+1)+1)]
                        for mode in range(0, M_R + 1):
                            Hm = jv(mode, k * L) + 1j * yv(mode, k * L)  # Hankel function type 1(+) or 2(-) ?
                            contribution = coef[M_R+mode] * Hm * np.exp(1j * mode * alpha)
                            if j==body : 
                                Hm_i = jv(mode, k * L_i) + 1j * yv(mode, k * L_i) 
                                contribution += coef_iso[M_R+mode] * Hm_i * np.exp(1j * mode * alpha_i)
                            
                            if mode>0:
                                contribution+=coef[M_R-mode] * Hm * np.exp(-1j * mode * alpha)*(-1)**(mode)
                                if j==body :
                                    contribution+=coef_iso[M_R-mode] * Hm_i * np.exp(-1j * mode * alpha_i)*(-1)**(mode)
                            if Ne>0:
                                for l in range(1, Ne+1) :
                                    # Km=kv(mode, kll[l-1]*L)
                                    if j==body :
                                        Km_i=kv(mode, kll[l-1]*L_i)
                                        coef_E_iso=(AR_iso[ind][dof][(2*M_R+1)*l:(2*M_R+1)*(l+1)])
                                        contribution += coef_E_iso[M_R+mode] * Km_i * np.cos(kll[l-1]*depth) * np.exp(1j * mode * alpha_i)
                                    # coef_E=(aR[ind][dof][body][(2*M_R+1)*(j*(Ne+1)+l):(2*M_R+1)*(j*(Ne+1)+l+1)])
                                    # contribution += coef_E[M_R+mode] * Km * np.cos(kll[l-1]*depth) * np.exp(1j * mode * alpha)
                                    if mode>0:
                                        # contribution+=coef_E[M_R-mode] * Km * np.cos(kll[l-1]*depth) * np.exp(-1j * mode * alpha)*(-1)**(mode)
                                        if j==body :
                                            contribution+=coef_E_iso[M_R-mode] * Km_i * np.cos(kll[l-1]*depth) * np.exp(-1j * mode * alpha_i)*(-1)**(mode)
                            
                            ETA_R[ind, index] += contribution   # for each dof and each body
                            phiR[ind] += contribution           # sum of every radiation problem
                                    
            # for diffraction problem
            for beta in range(Ndir):
                for j in range(Nb):  
                    dx = X - coord[j, 0]
                    dy = Y - coord[j, 1]
                    L = np.sqrt(dx**2 + dy**2)
                    alpha = np.arctan2(dy, dx)
                    coef=aS[ind][beta][(2*M_S+1)*j:(2*M_S+1)*(j+1)]
                    for mode in range(0, M_S + 1):
                        Hm = jv(mode, k * L) + 1j * yv(mode, k * L)
                        contribution = coef[M_S+mode] * Hm * np.exp(1j * mode * alpha)
                        if mode>0:
                            contribution+=coef[M_S-mode] * Hm * np.exp(-1j * mode * alpha)*(-1)**(mode)
                        # contribution = ((-1)**(-mode))*(g/w)*(1j*aS[ind][beta][idx_mode]) * Hm * np.exp(-1j * mode * alpha)
                        phiS[ind, beta] += contribution

            # normalization
            # *(1j*g/w) pour phi et *(-1j * w / g) pour ETA
            ETA_R[ind] *= (-1j*g/w) *(-1j*w/g) /(1j*w) #/2
            ETA_S[ind] = phiS[ind]*(1j*g/w)* (-1j * w / g)

            # calcul of the total potential
            for beta in range(Ndir):
                phiI = np.exp(1j * k * (X * np.cos(direction[beta]) + Y * np.sin(direction[beta])))
                ETA[ind, beta] = (phiS[ind, beta] + phiR[ind]+phiI) #* (1j * w / g)
            
        self.ETA=ETA
        self.ETA_R=ETA_R 
        self.ETA_S=ETA_S
    
   

    def Write(self, directory, Ne) : 

        AR=self.aR
        AS=self.aS
        w=2*np.pi/self.period
        if Ne>0 : 
            E=f"E{Ne}_"
        else : 
            E=""
        AR_file_path_abs = os.path.join(directory,  f"MultiBodyProblem_{E}RadiationCoefficients_abs.dat")
        AS_file_path_abs = os.path.join(directory,  f"MultiBodyProblem_{E}ScatteringCoefficients_abs.dat")
        AR_file_path_ph = os.path.join(directory,  f"MultiBodyProblem_{E}RadiationCoefficients_ph.dat")
        AS_file_path_ph = os.path.join(directory,  f"MultiBodyProblem_{E}ScatteringCoefficients_ph.dat")
        with open(AR_file_path_abs, "w") as AR_file:
                for i, period in enumerate(w):
                    for dof in range(len(AR[0])):  # Loop over the forces x bodies
                        for body in range(len(AR[0][0])):  # Loop over the forces x bodies
                            AR_file.write(f"{period:.4f}    ")
                            for mode in range(len(AR[i][0][0])):  # Loop over the forces x bodies
                                AR_file.write(f" {np.abs(AR[i][dof][body][mode]):.6e}  ")  # Absolute value of Fe
                            AR_file.write("\n")
        with open(AR_file_path_ph, "w") as AR_file:
                for i, period in enumerate(w):
                    for dof in range(len(AR[0])):  # Loop over the forces x bodies
                        for body in range(len(AR[0][0])):  # Loop over the forces x bodies
                            AR_file.write(f"{period:.4f}    ")
                            for mode in range(len(AR[i][0][0])):  # Loop over the forces x bodies
                                AR_file.write(f" {cmath.phase(AR[i][dof][body][mode]):.6e}  ")  # Absolute value of Fe
                            AR_file.write("\n")
        with open(AS_file_path_abs, "w") as AS_file:
                for i, period in enumerate(w):
                    for beta in range(len(AS[0])):  # Loop over the forces x bodies
                        AS_file.write(f"{period:.4f}    ")
                        for mode in range(len(AS[i][0])): # Loop over the forces x bodies
                            AS_file.write(f" {np.abs(AS[i][beta][mode]):.6e}  ")  # Absolute value of Fe
                        AS_file.write("\n")
        with open(AS_file_path_ph, "w") as AS_file:
                for i, period in enumerate(w):
                    for beta in range(len(AS[0])):  # Loop over the forces x bodies
                        AS_file.write(f"{period:.4f}    ")
                        for mode in range(len(AS[i][0])): # Loop over the forces x bodies
                            AS_file.write(f" {cmath.phase(AS[i][beta][mode]):.6e}  ")  # Absolute value of Fe
                        AS_file.write("\n")
      

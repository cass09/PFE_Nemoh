"""
"""

# Imports
from math import pi
import numpy as np
import os
from pickle import Pickler, Unpickler
from toolbox.CalaixSastre import WNumber, WNumber_E
from BEM import WAMIT, Nemoh
from semi_analytical.transfers import transfers
import cmath

# Settings
sep = os.path.join(' ' , ' ').split()[0]
TolCheck = 1e-3

class Body(object):

    def __init__(self, freqs, directions, depth, Nmodes_E, convention='W') :
        """
        """
        self.period = 2*pi/freqs
        self.dir = directions
        self.depth = depth
        self.Nmodes_E=Nmodes_E
        self.wnumber = WNumber(self.period, self.depth)
        self.wnumber_E = WNumber_E(self.period, self.depth, self.Nmodes_E)
        self.convention = convention

    def Transfers(self, dirSet, dirDy, dirFP, Evanescent, BEM ='N', Tol=1e-6, check=True):
        """
        """
        if BEM == 'W':
            # Read WAMIT output files
            (period, dirs, depth, NBODY, dof, dof_gener, dof_solid, dofi_solid) = WAMIT.ReadCFGandPOT(dirSet)
            (Fex, Madd, Crad) = WAMIT.ReadDynamics(dirDy, period, dirs, dof, self.convention)
            if dirFP:
                (PhiS, PhiR, r, t, z) = WAMIT.ReadFieldPoints(dirFP, period, dirs, depth, dof, dof_solid, self.convention)

        elif BEM == 'N':
            # Read Nemoh output files
            (dof, force, freq, dirs, FieldPoints, whichDOF, whichFORCE, depth) = Nemoh.ReadNemohcal(dirSet)
            period = 2*pi/freq
            (Fex, Madd, Crad) = Nemoh.ReadDynamics(dirDy, dof, force, len(dirs), self.convention)
            if dirFP:
                (PhiS, PhiR, r, t, z) = Nemoh.ReadFieldPoints(dirFP, dof, freq, len(dirs), FieldPoints, self.convention)

        # Assign new atributes
        self.Fex = Fex ; self.Madd = Madd ; self.Crad = Crad ; self.dof = dof
        if dirFP: self.radius = r ; self.fpazimuth = t ; self.fpdepth = z

        # Check point
        if check:
            requires = [abs(self.period-period) < TolCheck, abs(self.dir-dirs) < TolCheck, [abs(self.depth-depth) < TolCheck]]
            requires = [all(requires[ind]) for ind in range(len(requires))]
            if not all(requires):
                print('period {} , dirs {} , depth {}'.format(*requires))
                raise IOError('Wave data retrieved from output files is different than inputs')
        else:
            self.period, self.dir, self.depth = period, dirs, depth
            self.wnumber = WNumber(self.period, self.depth)
            self.wnumber_E = WNumber_E(self.period, self.depth, self.Nmodes_E)
        # Compute diffraction and force transfer matrices
        if dirFP:
            (self.D, self.G, self.AR, self.AS, self.order, self.truncorder) = transfers(self.depth, self.dir, self.period,
                                                                               (self.radius, self.fpazimuth, self.fpdepth),
                                                                               PhiS, PhiR, self.Fex, Tol, self.convention, 
                                                                               Evanescent, self.Nmodes_E)
        else:
            self.D = self.G = self.AR =self.AS = self.order = self.truncorder = self.radius = self.fpazimuth = self.fpdepth = dirFP

    def PickTransfers(self, picklefile, save=False, check=True) :
        """
        """
        # Generate a pickler
        if save :
            with open(picklefile , 'wb') as fid:
                iB = Pickler(fid).dump(self)
        else :
            with open(picklefile , 'rb') as fid:
                iB = Unpickler(fid).load()
            (self.Fex, self.Madd, self.Crad, self.radius, self.fpazimuth, self.fpdepth,
             depth, dirs, period, self.wnumber, self.dof, self.order, self.truncorder,
             self.G, self.D, self.AR, self.AS) = \
             (iB.Fex, iB.Madd, iB.Crad, iB.radius, iB.fpazimuth, iB.fpdepth,
             iB.depth, iB.dir, iB.period, iB.wnumber, iB.dof, iB.order, iB.truncorder,
             iB.G, iB.D, iB.AR, iB.AS)

            # Check point
            if check:
                requires = [abs(self.period-period) < TolCheck, abs(self.dir-dirs) < TolCheck, [abs(self.depth-depth) < TolCheck]]
                requires = [all(requires[ind]) for ind in range(len(requires))]
                if not all(requires):
                    print('period {} , dirs {} , depth {}'.format(*requires))
                    raise IOError('Wave data retrieved from output files is different than inputs')
            else: self.period, self.dir, self.depth = period, dirs, depth
        return iB

    def reduce_num_periods(self, per_red):
        """
        """
        cond = np.array([np.arange(len(self.period))[self.period == per][0] for per in per_red], int)
        self.period = self.period[cond]
        self.wnumber = self.wnumber[cond]
        self.wnumber_E = self.wnumber_E[cond]
        self.Fex = self.Fex[cond]
        self.Madd = self.Madd[cond]
        self.Crad = self.Crad[cond]
        self.D = self.D[cond]
        self.G = self.G[cond]
        self.AR = self.AR[cond]
        self.AS = self.AS[cond]
        self.truncorder = self.truncorder[cond]

    def Write(self, directory, Evanescent) : 

        D=self.D
        G=self.G
        AR=self.AR
        AS=self.AS
        w=2*np.pi/self.period
        if Evanescent : 
            E=f"E{self.Nmodes_E}_"
        else : 
            E=""
        D_file_path_abs = os.path.join(directory,  f"OneBodyProblem_{E}DiffractionMatrix_abs.dat")
        D_file_path_ph = os.path.join(directory,  f"OneBodyProblem_{E}DiffractionMatrix_ph.dat")
        D_file_path_Re = os.path.join(directory,  f"OneBodyProblem_{E}DiffractionMatrix_Re.dat")
        D_file_path_Imag = os.path.join(directory,  f"OneBodyProblem_{E}DiffractionMatrix_Imag.dat")
        G_file_path_abs = os.path.join(directory,  f"OneBodyProblem_{E}ForceTransferMatrix_abs.dat")
        G_file_path_ph = os.path.join(directory,  f"OneBodyProblem_{E}ForceTransferMatrix_ph.dat")
        AR_file_path_abs = os.path.join(directory,  f"OneBodyProblem_{E}RadiationCoefficients_abs.dat")
        AS_file_path_abs = os.path.join(directory,  f"OneBodyProblem_{E}ScatteringCoefficients_abs.dat")
        AR_file_path_ph = os.path.join(directory,  f"OneBodyProblem_{E}RadiationCoefficients_ph.dat")
        AR_file_path_Re = os.path.join(directory,  f"OneBodyProblem_{E}RadiationCoefficients_Re.dat")
        AR_file_path_Imag = os.path.join(directory,  f"OneBodyProblem_{E}RadiationCoefficients_Imag.dat")
        AS_file_path_ph = os.path.join(directory,  f"OneBodyProblem_{E}ScatteringCoefficients_ph.dat")
        with open(D_file_path_abs, "w") as D_file:
                for i, period in enumerate(w):
                    for j in range(len(D[0, :, 0])):  # Loop over the forces x bodies
                        D_file.write(f"{period:.4f}    ") 
                        for k in range(len(D[0, 0, :])):  # Loop over the forces x bodies
                            D_file.write(f" {np.abs(D[i, j, k]):.6e}  ")  # Absolute value of Fe
                        D_file.write("\n")
        with open(D_file_path_ph, "w") as D_file:
                for i, period in enumerate(w):
                    for j in range(len(D[0, :, 0])):  # Loop over the forces x bodies
                        D_file.write(f"{period:.4f}    ") 
                        for k in range(len(D[0, 0, :])):  # Loop over the forces x bodies
                            D_file.write(f" {cmath.phase(D[i, j, k]):.6e}  ")  # Absolute value of Fe
                        D_file.write("\n")
        with open(D_file_path_Re, "w") as D_file:
                for i, period in enumerate(w):
                    for j in range(len(D[0, :, 0])):  # Loop over the forces x bodies
                        D_file.write(f"{period:.4f}    ") 
                        for k in range(len(D[0, 0, :])):  # Loop over the forces x bodies
                            D_file.write(f" {(D[i, j, k]).real:.6e}  ")  # Absolute value of Fe
                        D_file.write("\n")
        with open(D_file_path_Imag, "w") as D_file:
                for i, period in enumerate(w):
                    for j in range(len(D[0, :, 0])):  # Loop over the forces x bodies
                        D_file.write(f"{period:.4f}    ") 
                        for k in range(len(D[0, 0, :])):  # Loop over the forces x bodies
                            D_file.write(f" {(D[i, j, k]).imag:.6e}  ")  # Absolute value of Fe
                        D_file.write("\n")
        with open(G_file_path_abs, "w") as G_file:
                for i, period in enumerate(w):
                    for j in range(len(G[0, :, 0])):  # Loop over the forces x bodies
                        G_file.write(f"{period:.4f}    ")
                        for k in range(len(G[0, 0, :])):  # Loop over the forces x bodies
                            G_file.write(f" {np.abs(G[i, j, k]):.6e}  ")  # Absolute value of Fe
                        G_file.write("\n")
        with open(G_file_path_ph, "w") as G_file:
                for i, period in enumerate(w):
                    for j in range(len(G[0, :, 0])):  # Loop over the forces x bodies
                        G_file.write(f"{period:.4f}    ")
                        for k in range(len(G[0, 0, :])):  # Loop over the forces x bodies
                            G_file.write(f" {cmath.phase(G[i, j, k]):.6e}  ")  # Absolute value of Fe
                        G_file.write("\n")
        with open(AR_file_path_abs, "w") as AR_file:
                for i, period in enumerate(w):
                    for k in range(len(AR[0, 0, :])):  # Loop over the forces x bodies
                        AR_file.write(f"{period:.4f}    ")
                        for j in range(len(AR[0, :, 0])):  # Loop over the forces x bodies
                            AR_file.write(f" {np.abs(AR[i, j, k]):.6e}  ")  # Absolute value of Fe
                        AR_file.write("\n")
        with open(AR_file_path_ph, "w") as AR_file:
                for i, period in enumerate(w):
                    for k in range(len(AR[0, 0, :])):  # Loop over the forces x bodies
                        AR_file.write(f"{period:.4f}    ")
                        for j in range(len(AR[0, :, 0])):  # Loop over the forces x bodies
                            AR_file.write(f" {cmath.phase(AR[i, j, k]):.6e}  ")  # Absolute value of Fe
                        AR_file.write("\n")
        with open(AR_file_path_Re, "w") as AR_file:
                for i, period in enumerate(w):
                    for k in range(len(AR[0, 0, :])):  # Loop over the forces x bodies
                        AR_file.write(f"{period:.4f}    ")
                        for j in range(len(AR[0, :, 0])):  # Loop over the forces x bodies
                            AR_file.write(f" {(AR[i, j, k]).real:.6e}  ")  # Absolute value of Fe
                        AR_file.write("\n")
        with open(AR_file_path_Imag, "w") as AR_file:
                for i, period in enumerate(w):
                    for k in range(len(AR[0, 0, :])):  # Loop over the forces x bodies
                        AR_file.write(f"{period:.4f}    ")
                        for j in range(len(AR[0, :, 0])):  # Loop over the forces x bodies
                            AR_file.write(f" {(AR[i, j, k]).imag:.6e}  ")  # Absolute value of Fe
                        AR_file.write("\n")
        with open(AS_file_path_abs, "w") as AS_file:
                for i, period in enumerate(w):
                    for k in range(len(AS[0, 0, :])):  # Loop over the forces x bodies
                        AS_file.write(f"{period:.4f}    ")
                        for j in range(len(AS[0, :, 0])):  # Loop over the forces x bodies
                            AS_file.write(f" {np.abs(AS[i, j, k]):.6e}  ")  # Absolute value of Fe
                        AS_file.write("\n")
        with open(AS_file_path_ph, "w") as AS_file:
                for i, period in enumerate(w):
                    for k in range(len(AS[0, 0, :])):  # Loop over the forces x bodies
                        AS_file.write(f"{period:.4f}    ")
                        for j in range(len(AS[0, :, 0])):  # Loop over the forces x bodies
                            AS_file.write(f" {cmath.phase(AS[i, j, k]):.6e}  ")  # Absolute value of Fe
                        AS_file.write("\n")
"""
"""

# Imports
from math import pi
import numpy as np
import os
from pickle import Pickler, Unpickler
from capytaine.interaction_theory.transfers_sources import transfers_sources
from capytaine.interaction_theory.transfers import transfers
import cmath

# Settings
sep = os.path.join(' ' , ' ').split()[0]
TolCheck = 1e-3

class Body(object):

    def __init__(self, freqs, directions, depth, Nmodes_E, wavenumber, Fex, MethodSource) :
        """
        """
        self.period = 2*pi/freqs
        self.dir = directions
        self.depth = depth
        self.Nmodes_E=Nmodes_E
        self.wnumber = wavenumber
        self.Fex=Fex
        self.method=MethodSource

    def Transfers(self, ParamS, ParamR, mesh, Tol=1e-6):
        if self.method:
            (self.D, self.G, self.AR, self.AS, self.order, self.truncorder) = transfers_sources(self.depth, self.dir, 
                self.period, self.wnumber, mesh, ParamS, ParamR, self.Fex, Tol, self.Nmodes_E)
        else : 
            (self.D, self.G, self.AR, self.AS, self.order, self.truncorder) = transfers(self.depth, self.dir, 
                self.period, self.wnumber, mesh, ParamS, ParamR, self.Fex, Tol, self.Nmodes_E)
       

    def Write(self, directory, L, method) : 

        D=self.D
        G=self.G
        AR=self.AR
        AS=self.AS
        w=2*np.pi/self.period
        if L>0 : 
            E=f"E{self.Nmodes_E}_"
        else : 
            E=""
        if method==1:
            S="S_"
        else:
            S=""
        D_file_path_abs = os.path.join(directory,  f"{S}OneBodyProblem_{E}DiffractionMatrix_abs.dat")
        D_file_path_ph = os.path.join(directory,  f"{S}OneBodyProblem_{E}DiffractionMatrix_ph.dat")
        D_file_path_Re = os.path.join(directory,  f"{S}OneBodyProblem_{E}DiffractionMatrix_Re.dat")
        D_file_path_Imag = os.path.join(directory,  f"{S}OneBodyProblem_{E}DiffractionMatrix_Imag.dat")
        G_file_path_Re = os.path.join(directory,  f"{S}OneBodyProblem_{E}ForceTransferMatrix_Re.dat")
        G_file_path_Imag = os.path.join(directory,  f"{S}OneBodyProblem_{E}ForceTransferMatrix_Imag.dat")
        AR_file_path_abs = os.path.join(directory,  f"{S}OneBodyProblem_{E}RadiationCoefficients_abs.dat")
        AS_file_path_abs = os.path.join(directory,  f"{S}OneBodyProblem_{E}ScatteringCoefficients_abs.dat")
        AR_file_path_ph = os.path.join(directory,  f"{S}OneBodyProblem_{E}RadiationCoefficients_ph.dat")
        AR_file_path_Re = os.path.join(directory,  f"{S}OneBodyProblem_{E}RadiationCoefficients_Re.dat")
        AR_file_path_Imag = os.path.join(directory,  f"{S}OneBodyProblem_{E}RadiationCoefficients_Imag.dat")
        AS_file_path_ph = os.path.join(directory,  f"{S}OneBodyProblem_{E}ScatteringCoefficients_ph.dat")
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
        with open(G_file_path_Re, "w") as G_file:
                for i, period in enumerate(w):
                    for j in range(len(G[0, :, 0])):  # Loop over the forces x bodies
                        G_file.write(f"{period:.4f}    ")
                        for k in range(len(G[0, 0, :])):  # Loop over the forces x bodies
                            G_file.write(f" {(G[i, j, k]).real:.6e}  ")  # Absolute value of Fe
                        G_file.write("\n")
        with open(G_file_path_Imag, "w") as G_file:
                for i, period in enumerate(w):
                    for j in range(len(G[0, :, 0])):  # Loop over the forces x bodies
                        G_file.write(f"{period:.4f}    ")
                        for k in range(len(G[0, 0, :])):  # Loop over the forces x bodies
                            G_file.write(f" {(G[i, j, k]).imag:.6e}  ")  # Absolute value of Fe
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
        with open(os.path.join(directory,f"Data_waves.dat"), "w") as f:
            f.write(f"w (rad/s) - k (m^-1) - L or lambda (m) - T (s) - f (Hz)\n  ")
            for i, value in enumerate(w) : 
                f.write(f"{value:8.4f}  {self.wnumber[i]:8.4f}    {2*pi/self.wnumber[i]:8.4f}     {2*np.pi/value:8.4f}  {value/(2*np.pi):8.4f}\n  ")
        # with open(os.path.join(directory,f"Data_Ewaves.dat"), "w") as f:
        #     f.write(f"w (rad/s) - kl (m^-1) with L={self.Nmodes_E} \n")
        #     for i, value in enumerate(w) : 
        #         f.write(f"{value:8.4f}")
        #         for e in range(self.Nmodes_E) :
        #             f.write(f" {self.wnumber_E[i, e]:8.4f} ")
        #         f.write("\n ")

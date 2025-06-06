import numpy as np
import cmath
import os
import h5py

def WriteData(WECArr, directionMB, farm, distance, results, results_h5, Motion):
    
    w = 2*np.pi/WECArr.period
    FR = 1j*w.T*WECArr.Madd.T + WECArr.Crad.T
    if len(directionMB)==1 and directionMB[0]!=0.00 :
        beta_value=f"beta_{directionMB[0]}_"
    else :
        beta_value=""
    if farm['Configuration'] :
        type=farm["Type"] + "_"
    else :
        type=""
    if farm['Ne_modes']>0 : 
        test= f"E{farm['Ne_modes']}_"
    else : 
        test=""

    if farm["Output_type"]=="DAT" : 
        fe_abs_file_path = os.path.join(results,  f"Global_{test}Fe_abs_{type}{beta_value}d{distance:.2f}.dat")
        fe_phase_file_path = os.path.join(results,  f"Global_{test}Fe_phase_{type}{beta_value}d{distance:.2f}.dat")
        madd_file_path = os.path.join(results, f"Global_{test}Madd_{type}{beta_value}d{distance:.2f}.dat")
        crad_file_path = os.path.join(results,  f"Global_{test}Crad_{type}{beta_value}d{distance:.2f}.dat")

        
        #################################################################
        #################### Excitation Forces ##########################
        #################################################################
        with open(fe_abs_file_path, "w") as fe_file, open(fe_phase_file_path, "w") as fe_phase_file, open(madd_file_path, "w") as madd_file, open(crad_file_path, "w") as crad_file:
            fe_file.write("Period   " + " Wave direction    ")
            fe_file.write(" Fe dof = body x force")
            fe_file.write("\n")
            for k, beta in enumerate(directionMB):
                for i, period in enumerate(w):
                    fe_file.write(f"{beta:.{4}e}  ")
                    fe_file.write(f"{period:.{4}e}    ")  # Start with the period
                    for j in range(len(WECArr.Fex[0, 0, :])):  # Loop over the forces x bodies
                        fe_file.write(f" {np.abs(WECArr.Fex[i, k, j]):.6e}  ")   
                    fe_file.write("\n")
                fe_file.write("\n")

            fe_phase_file.write(" Wave direction   " + " Period    ")
            fe_phase_file.write(" Fe dof = body x force")
            fe_phase_file.write("\n")
            for k, beta in enumerate(directionMB):
                for i, period in enumerate(w):
                    fe_phase_file.write(f"{beta:.4f}  ")
                    fe_phase_file.write(f"{period:.4f}    ")  # Start with the period
                    for j in range(len(WECArr.Fex[0, 0, :])):  # Loop over the forces x bodies
                        fe_phase_file.write(f" {cmath.phase(WECArr.Fex[i, k, j]):.6e}  ")  # Phase of Fe
                    fe_phase_file.write("\n")
                fe_phase_file.write("\n")

        #################################################################
        ################## Added Mass and Damping #######################
        #################################################################
            N = len(WECArr.Madd[0, :, 0])  # N est la dimension de la matrice Madd
            madd_file.write("Period " + "Madd_{ij}" + "\n")
            for i, period in enumerate(w):
                for j in range(N):
                    madd_file.write(f"{period:.4f}  ")  # Start with the period
                    for k in range(N):
                        madd_file.write(f" {WECArr.Madd[i, j, k]:.6e}  ")  # Madd values
                    madd_file.write("\n")

            crad_file.write("Period " + "Crad_{ij}" + "\n")
            for i, period in enumerate(w):
                for j in range(N):
                    crad_file.write(f"{period:.4f} ")  # Start with the period
                    for k in range(N):
                        crad_file.write(f" {WECArr.Crad[i, j, k]:.6e}  ")  # Crad values
                    crad_file.write("\n")

        #################################################################
        ########################### RAO #################################
        #################################################################
        if farm['RAO'] :
            RAO_file_path = os.path.join(Motion,  f"Global_{test}RAO_{type}{beta_value}d{distance:.2f}.dat")
            with open(RAO_file_path, "w") as RAO_file:

                RAO_file.write(" Frecency   "+"|X| (m/m)" +"|Y| (m/m)"+ "|Z| (m/m)" +"|phi| (deg)" +"|theta| (deg)" +"|psi| (deg)" +"ang(x) (deg)"+ "ang(y) (deg)"+ "ang(z) (deg)"+ "ang(phi) (deg)"+ "ang(theta) (deg)" +"ang(psi) (deg)")
                RAO_file.write("\n")
                for k, beta in enumerate(directionMB):
                    RAO_file.write(f"beta={beta:.4f} \n ")
                    for i, period in enumerate(w):
                        RAO_file.write(f"{period:.4f}    ")  # Start with the period
                        for j in range(len(WECArr.RAO[0, 0, :])):  # Loop over the forces x bodies
                            RAO_file.write(f" {np.abs(WECArr.RAO[i, k, j]):.6e}  ")   
                        for j in range(len(WECArr.RAO[0, 0, :])):  # Loop over the forces x bodies
                            RAO_file.write(f" {cmath.phase(WECArr.RAO[i, k, j]):.6e}  ")   
                        RAO_file.write("\n")

        #################################################################
        ######################### Kochin ################################
        #################################################################
        if farm['Kochin']['number']>0 :
            thetaK = np.linspace(farm['Kochin']['min'],
                        farm['Kochin']['max'],
                        farm['Kochin']['number'])
            if farm['Kochin']['format']=='DEG' :
                thetaK=thetaK*np.pi/180
            ind_pb=1
            for i, period in enumerate(w): 
                for dir in range(len(WECArr.KochinS[0, :, 0])):
                    Kochin_file_path = os.path.join(Motion,  f"Global_{test}Kochin_{type}{beta_value}d{distance:.2f}_pb{ind_pb}.dat")
                    with open(Kochin_file_path, "w") as K_file:
                        for theta in range(len(WECArr.KochinS[0, 0, :])):  
                            K_file.write(f"     {thetaK[theta]:.6e}  ")      
                            K_file.write(f" {np.abs(WECArr.KochinS[i, dir, theta]):.6e}  ")      
                            K_file.write(f" {cmath.phase(WECArr.KochinS[i, dir, theta]):.6e}  ")   
                            K_file.write(f" {(WECArr.KochinS[i, dir, theta]).real:.6e}  ")   
                            K_file.write(f" {(WECArr.KochinS[i, dir, theta]).imag:.6e}  ")   
                            K_file.write("\n")
                    ind_pb=ind_pb+1
                for dof in range(len(WECArr.KochinR[0, :, 0])): 
                    Kochin_file_path = os.path.join(Motion,  f"Global_{test}Kochin_{type}{beta_value}d{distance:.2f}_pb{ind_pb}.dat")
                    with open(Kochin_file_path, "w") as K_file:
                        for theta in range(len(WECArr.KochinR[0, 0, :])): 
                            K_file.write(f"     {thetaK[theta]:.6e}  ")   
                            K_file.write(f" {np.abs(WECArr.KochinR[i, dof, theta]):.6e}  ")      
                            K_file.write(f" {cmath.phase(WECArr.KochinR[i, dof, theta]):.6e}  ")   
                            K_file.write(f" {(WECArr.KochinR[i, dof, theta]).real:.6e}  ")   
                            K_file.write(f" {(WECArr.KochinR[i, dof, theta]).imag:.6e}  ")   
                            K_file.write("\n")
                    ind_pb=ind_pb+1

        #################################################################
        ###################### Free surface #############################
        #################################################################
        if farm['Free_surface']['Nx'] > 0 :
            ind_pb=1
            Nx=farm['Free_surface']['Nx']
            Ny=farm['Free_surface']['Ny']
            Lx=farm['Free_surface']['Ly']
            Ly=farm['Free_surface']['Ly']
            x = np.linspace(-Lx/2, Lx/2, Nx)
            y = np.linspace(-Ly/2, Ly/2, Ny)
            for ind, period in enumerate(w): 
                for dir in range(len(WECArr.ETA_S[0, :, 0, 0])):
                    FS_file_path = os.path.join(Motion,  f"Global_{test}FS_{type}{beta_value}d{distance:.2f}_pb{ind_pb:05d}.dat")
                    with open(FS_file_path, "w") as FS_file:
                        for i in range(Nx):  
                            for j in range(Ny):  
                                FS_file.write(f"     {x[i]:.6e}  ")   
                                FS_file.write(f"     {y[j]:.6e}  ")       
                                FS_file.write(f" {np.abs(WECArr.ETA_S[ind, dir, i,j]):.6e}  ")      
                                FS_file.write(f" {cmath.phase(WECArr.ETA_S[ind, dir, i,j]):.6e}  ")   
                                FS_file.write(f" {(WECArr.ETA_S[ind, dir, i,j]).real:.6e}  ")   
                                FS_file.write(f" {(WECArr.ETA_S[ind, dir, i,j]).imag:.6e}  ")   
                                FS_file.write("\n")
                    ind_pb=ind_pb+1
                for dof in range(len(WECArr.ETA_R[0, :, 0, 0])): 
                    FS_file_path = os.path.join(Motion,  f"Global_{test}FS_{type}{beta_value}d{distance:.2f}_pb{ind_pb:05d}.dat")
                    with open(FS_file_path, "w") as FS_file:
                        for i in range(Nx):  
                            for j in range(Ny):  
                                FS_file.write(f"     {x[i]:.6e}  ")   
                                FS_file.write(f"     {y[j]:.6e}  ")       
                                FS_file.write(f" {np.abs(WECArr.ETA_R[ind, dir, i,j]):.6e}  ")      
                                FS_file.write(f" {cmath.phase(WECArr.ETA_R[ind, dir, i,j]):.6e}  ")   
                                FS_file.write(f" {(WECArr.ETA_R[ind, dir, i,j]).real:.6e}  ")   
                                FS_file.write(f" {(WECArr.ETA_R[ind, dir, i,j]).imag:.6e}  ")   
                                FS_file.write("\n")
                    ind_pb=ind_pb+1
                    
            ind_pb=1  
            for ind, period in enumerate(w): 
                for dir in range(len(WECArr.ETA[0, :, 0, 0])):
                    ETA_file_path = os.path.join(Motion,  f"Global_{test}ETA_{type}{beta_value}d{distance:.2f}_pb{ind_pb:05d}.dat")
                    with open(ETA_file_path, "w") as FS_file:
                        FS_file.write(f"w={period:.4f}  ")  
                        FS_file.write(f"beta={directionMB[dir]:.4f}  ")  
                        FS_file.write("\n")  
                        for i in range(Nx):  
                            for j in range(Ny):  
                                FS_file.write(f"     {x[i]:.6e}  ")   
                                FS_file.write(f"     {y[j]:.6e}  ")       
                                FS_file.write(f" {np.abs(WECArr.ETA[ind, dir, i,j]):.6e}  ")      
                                FS_file.write(f" {cmath.phase(WECArr.ETA[ind, dir, i,j]):.6e}  ")   
                                FS_file.write(f" {(WECArr.ETA[ind, dir, i,j]).real:.6e}  ")   
                                FS_file.write(f" {(WECArr.ETA[ind, dir, i,j]).imag:.6e}  ")   
                                FS_file.write("\n")
                    ind_pb=ind_pb+1

    else :
        h5f = h5py.File(results_h5, 'w')
        h5f.create_dataset(f'radiation_{type}_d{distance}', data=FR)
        h5f.create_dataset(f'excitation_{type}_d{distance}', data=WECArr.Fex)
        h5f.create_dataset('frequency', data=WECArr.period)
        h5f.close()

    return farm["Output_type"]
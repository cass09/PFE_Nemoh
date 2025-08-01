import os
import numpy as np
import cmath

def WriteCapytaineData(dataset, rao_dataset, results_folder, dofs=None):
    """
    Exporte les résultats de Capytaine (Fex, Madd, Crad, RAO) dans des fichiers .dat lisibles.
    """

    # Crée le dossier si besoin
    os.makedirs(results_folder, exist_ok=True)

    omegas = dataset.coords['omega'].values
    periods = 2 * np.pi / omegas
    wave_directions = dataset.coords['wave_direction'].values
    radiating_dofs = dataset.coords['radiating_dof'].values
    influenced_dof = dataset.coords['influenced_dof'].values

    if dofs is not None:
        dof_mask = [dof in dofs for dof in radiating_dofs]
    else:
        dof_mask = [True] * len(radiating_dofs)

    #################################################################
    #################### Excitation Forces ##########################
    #################################################################
    fe_abs_file = os.path.join(results_folder, "Capytaine_Fe_abs.dat")
    fe_phase_file = os.path.join(results_folder, "Capytaine_Fe_phase.dat")
    dataset['excitation_force'] = dataset['Froude_Krylov_force'] + dataset['diffraction_force']

    with open(fe_abs_file, "w") as fe_file, open(fe_phase_file, "w") as fe_phase:
        fe_file.write("Frequency     direction     " + "   ".join([d for i, d in enumerate(radiating_dofs) if dof_mask[i]]) + "\n")
        fe_phase.write("Frequency     direction     " + "   ".join([d for i, d in enumerate(radiating_dofs) if dof_mask[i]]) + "\n")

        for beta in wave_directions:
            for omega in omegas:
                fe_file.write(f"{beta:.6e}  {omega:.6e}  ")
                fe_phase.write(f"{beta:.6e}  {omega:.6e}  ")
                for i, dof in enumerate(influenced_dof):
                    if not dof_mask[i]:
                        continue
                    Fex = dataset["excitation_force"].sel(omega=omega, wave_direction=beta, influenced_dof=dof).values
                    fe_file.write(f"{np.abs(Fex):.6e}  ")
                    fe_phase.write(f"{np.angle(Fex):.6e}  ")
                fe_file.write("\n")
                fe_phase.write("\n")

    #################################################################
    ##################### Added Mass & Damping ######################
    #################################################################
    madd_file = os.path.join(results_folder, "Capytaine_Madd.dat")
    crad_file = os.path.join(results_folder, "Capytaine_Crad.dat")

    N = len(radiating_dofs)

    with open(madd_file, "w") as madd_f, open(crad_file, "w") as crad_f:
        madd_f.write(f"Frequency  {N} dof\n")
        crad_f.write(f"Frequency  {N} dof\n")

        for i_omega, omega in enumerate(omegas):
            for i in range(N):
                madd_f.write(f"{omega:.4f}  ")
                crad_f.write(f"{omega:.4f}  ")
                for j in range(N):
                    added_mass = dataset["added_mass"].sel(omega=omega,
                                                           radiating_dof=radiating_dofs[i],
                                                           influenced_dof=radiating_dofs[j]).values
                    damping = dataset["radiation_damping"].sel(omega=omega,
                                                               radiating_dof=radiating_dofs[i],
                                                               influenced_dof=radiating_dofs[j]).values
                    madd_f.write(f"{added_mass:.6e}  ")
                    crad_f.write(f"{damping:.6e}  ")
                madd_f.write("\n")
                crad_f.write("\n")

    #################################################################
    ########################### RAO #################################
    #################################################################
    rao_file = os.path.join(results_folder, "Capytaine_RAO.dat")
    with open(rao_file, "w") as rao_f:
        dofs_list = list(radiating_dofs)
        rao_f.write("Frequency   " + " ".join([f"|{dof}|(m/m)" for dof in dofs_list]) + "  " +
                    " ".join([f"ang({dof})(deg)" for dof in dofs_list]) + "\n")

        for beta in wave_directions:
            rao_f.write(f"beta = {beta:.4f}\n")
            for i_omega, omega in enumerate(omegas):
                rao_f.write(f"{omega:.4f}  ")
                for i, dof in enumerate(radiating_dofs):
                    if not dof_mask[i]:
                        continue
                    val = rao_dataset.sel(omega=omega, wave_direction=beta, radiating_dof=dof).values
                    rao_f.write(f"{np.abs(val):.6e}  ")
                for i, dof in enumerate(radiating_dofs):
                    if not dof_mask[i]:
                        continue
                    val = rao_dataset.sel(omega=omega, wave_direction=beta, radiating_dof=dof).values
                    rao_f.write(f"{np.degrees(np.angle(val)):.6e}  ")
                rao_f.write("\n")
        

def WriteCapytaineDataFromResults(results_data, results_folder, dofs=None):
    """
    Exporte les résultats (résolus pb par pb) dans des fichiers .dat :
    - Forces d'excitation
    - Masse ajoutée
    - Amortissement de radiation
    - RAOs
    """


    os.makedirs(results_folder, exist_ok=True)

    # Préparation des fichiers
    fe_abs_file = os.path.join(results_folder, "Capytaine_Fe_abs.dat")
    fe_phase_file = os.path.join(results_folder, "Capytaine_Fe_phase.dat")
    madd_file = os.path.join(results_folder, "Capytaine_Madd.dat")
    crad_file = os.path.join(results_folder, "Capytaine_Crad.dat")
    rao_file = os.path.join(results_folder, "Capytaine_RAO.dat")

    # Trouver les dofs
    all_dofs = sorted({res["dof"] for res in results_data if res["dof"] is not None})
    if dofs is None:
        dofs = all_dofs
    dof_mask = [dof in dofs for dof in all_dofs]

    # ---------- Excitation Forces ----------
    with open(fe_abs_file, "w") as fe_file, open(fe_phase_file, "w") as fe_phase:
        header = "direction     Frequency     " + "   ".join(dofs) + "\n"
        fe_file.write(header)
        fe_phase.write(header)

        for res in results_data:
            if res["type"] != "diffraction":
                continue
            omega = res["omega"]
            beta = res["beta"]
            fe_file.write(f"{beta:.6e}  {omega:.6e}  ")
            fe_phase.write(f"{beta:.6e}  {omega:.6e}  ")

            for dof in dofs:
                Fex = res["excitation_force"].get(dof, 0.0)
                fe_file.write(f"{np.abs(Fex):.6e}  ")
                fe_phase.write(f"{np.angle(Fex):.6e}  ")
            fe_file.write("\n")
            fe_phase.write("\n")

    # ---------- Added Mass & Damping ----------
    with open(madd_file, "w") as madd_f, open(crad_file, "w") as crad_f:
        madd_f.write(f"Frequency  {len(dofs)} dof\n")
        crad_f.write(f"Frequency  {len(dofs)} dof\n")

        unique_omegas = sorted(set(res["omega"] for res in results_data))
        for omega in unique_omegas:
            for dof_i in dofs:
                madd_f.write(f"{omega:.4f}  ")
                crad_f.write(f"{omega:.4f}  ")
                for dof_j in dofs:
                    val = next((res for res in results_data if
                                res["type"] == "radiation" and
                                res["omega"] == omega and
                                res["dof"] == dof_i), None)
                    if val is not None:
                        madd = val["added_mass"].get(dof_j, 0.0)
                        damp = val["damping"].get(dof_j, 0.0)
                    else:
                        madd = damp = 0.0
                    madd_f.write(f"{madd:.6e}  ")
                    crad_f.write(f"{damp:.6e}  ")
                madd_f.write("\n")
                crad_f.write("\n")


def WriteCapytaineRAO(rao, dofs, output_folder):
    """
    Écrit les RAOs calculés depuis `compute_rao_from_results()` dans un fichier texte lisible.

    Parameters
    ----------
    rao : dict
        Résultat retourné par `compute_rao_from_results()`, c’est-à-dire rao[omega][beta] = array (complex).
    dofs : list of str
        Liste des DOFs dans l’ordre utilisé pour construire les vecteurs RAO.
    output_file : str
        Chemin du fichier de sortie.
    """
    omegas = sorted(rao.keys())
    betas = sorted({beta for w in rao for beta in rao[w]})
    output_file = os.path.join(output_folder, "Capytaine_RAO.dat")
    with open(output_file, "w") as rao_f:
        # Écrire l'en-tête
        rao_f.write("Frequency   " + " ".join([f"|{dof}|(m/m)" for dof in dofs]) + "  " +
                    " ".join([f"ang({dof})(deg)" for dof in dofs]) + "\n")

        for beta in betas:
            rao_f.write(f"\nbeta = {beta:.4f}\n")
            for w in omegas:
                if beta not in rao[w]:
                    continue  # Pas de donnée pour ce couple (w, beta)
                rao_vec = rao[w][beta]
                rao_f.write(f"{w:.4f}  ")
                # Modules
                for val in rao_vec:
                    rao_f.write(f"{np.abs(val):.6e}  ")
                # Phases (en degrés)
                for val in rao_vec:
                    rao_f.write(f"{np.degrees(np.angle(val)):.6e}  ")
                rao_f.write("\n")



def WriteSourcesByProblem(results_data, body, results_folder="./sources"):
    os.makedirs(results_folder, exist_ok=True)

    collocation_points = body.mesh.faces_centers  # ou vertices, selon la formulation

    for i_pb, result in enumerate(results_data):
        sources_sigma = result.get("sources", None)
        if sources_sigma is None:
            print(f"[WARNING] No sources found for problem {i_pb+1}")
            continue
        
        if len(sources_sigma) != len(collocation_points):
            print(f"[ERROR] Mismatch between number of sources and collocation points in problem {i_pb+1}")
            continue
            
        filename = os.path.join(results_folder, f"sources.{i_pb+1:05d}.dat")
        with open(filename, "w") as f:
            for point, sigma in zip(collocation_points, sources_sigma):
                x, y, z = point
                re, im = np.real(sigma), np.imag(sigma)
                f.write(f"{re:.6e} {im:.6e}\n")
                # f.write(f"{x:.6e} {y:.6e} {z:.6e} {re:.6e} {im:.6e}\n")
                
def WriteCylsurfaceByProblem(results_data, OCmesh, results_folder="./cylsurface"):
    os.makedirs(results_folder, exist_ok=True)

    # collocation_points = body.mesh.faces_centers  # ou vertices, selon la formulation

    for i_pb, result in enumerate(results_data):
        cylsurface = result.get("cylsurface", None)
        if cylsurface is None:
            print(f"[WARNING] No cylsurface found for problem {i_pb+1}")
            continue
        
        if len(cylsurface) != len(OCmesh):
            print(f"[ERROR] Mismatch between number of sources and collocation points in problem {i_pb+1}")
            continue
        # Vérifie si c'est un problème de radiation
        typePB=result["type"]
        omega = result["omega"]

        filename = os.path.join(results_folder, f"cylsurface.{i_pb+1:05d}.dat")
        with open(filename, "w") as f:
            for point, sigma in zip(OCmesh, cylsurface):
                # if typePB=="radiation":
                #     sigma *= 1j / omega  # Modification pour les problèmes de radiation
                x, y, z = point
                re, im = np.real(sigma), np.imag(sigma)
                abs, ph = np.abs(sigma), np.angle(sigma)
                # f.write(f"{re:.6e} {im:.6e}\n")
                f.write(f"{x:.6e} {y:.6e} {z:.6e} {abs:.6e} {ph:.6e} {re:.6e} {im:.6e}\n")

def write_params_file(results_data, water_depth, output_file="params.dat"):
    with open(output_file, "w") as f:
        # En-tête
        f.write(f"depth={water_depth}\n")
        f.write("num\ttype\tomega\tbeta\tdof\tmodeM\tmodeL\n")
        for num, result in enumerate(results_data):
            type_ = result.get("type", "NA")
            
            omega = result.get("omega")
            omega_str = f"{omega:.6f}" if omega is not None else "NA"
            
            beta = result.get("beta")
            beta_str = f"{beta:.6f}" if beta is not None else "NA"
            
            dof = result.get("dof", "NA")
            modeM = result.get("modeM", "NA")
            modeL = result.get("modeL", "NA")
            
            # Écriture ligne
            f.write(f"{num+1}\t{type_}\t{omega_str}\t{beta_str}\t{dof}\t{modeM}\t{modeL}\n")



def extract_hydrodynamic_quantities_with_beta(results_data, dofs, unique_omegas, unique_betas, L):
    Nw = len(unique_omegas)
    Nbeta = len(unique_betas)
    Ndof = len(dofs)

    Madd = np.zeros((Nw, Ndof, Ndof), dtype=float)
    Crad = np.zeros((Nw, Ndof, Ndof), dtype=float)
    Fex = np.zeros((Nw, Nbeta*(L+1), Ndof), dtype=complex)

    for w_idx, omega in enumerate(unique_omegas):
        for i, dof_i in enumerate(dofs):
            # --- RADIATION ---
            val_rad = next((res for res in results_data if
                            res["type"] == "radiation" and
                            res["omega"] == omega and
                            res["dof"] == dof_i), None)
            if val_rad is not None:
                for j, dof_j in enumerate(dofs):
                    Madd[w_idx, i, j] = val_rad["added_mass"].get(dof_j, 0.0)
                    Crad[w_idx, i, j] = val_rad["damping"].get(dof_j, 0.0)

        # --- DIFFRACTION (Fex pour chaque beta) ---
        for l in range(L+1):  # mode vertical
            for b_idx, beta in enumerate(unique_betas):  # mode angulaire
                val_diff = next((res for res in results_data if
                                res["type"] == "diffraction" and
                                res["omega"] == omega and
                                res["beta"] == beta and
                                res.get("modeL", 0) == l), None)
                if val_diff is not None and isinstance(val_diff["excitation_force"], dict):
                    idx = l * Nbeta + b_idx
                    for i, dof in enumerate(dofs):
                        Fex[w_idx, idx, i] = val_diff["excitation_force"].get(dof, 0.0)

    return Madd, Crad, Fex




def extract_sources_by_problem(results_data, unique_omegas, unique_betas, dofs, L, Npanels):
    Nw = len(unique_omegas)
    Nb = len(unique_betas)
    Ndof = len(dofs)
    SourcesD = np.zeros((Nw, Nb*(L+1), Npanels), dtype=complex)
    SourcesR = np.zeros((Nw, Ndof, Npanels), dtype=complex)
    SourcesRb = np.zeros((Nw, Ndof, Npanels), dtype=complex)

    for w_idx, omega in enumerate(unique_omegas):
        # Diffraction
        for l in range(L + 1):                # modes verticaux
            for b_idx, beta in enumerate(unique_betas):  # modes angulaires m
                res_diff = next((res for res in results_data if
                                    res["type"] == "diffraction" and
                                    res["omega"] == omega and
                                    res["beta"] == beta and
                                    res.get("modeL", 0) == l), None)
                if res_diff is not None:
                    sigma = np.array(res_diff["sources"])
                    # Calcul de l’indice dans la deuxième dimension : décalage selon l et m
                    idx = l * Nb + b_idx
                    SourcesD[w_idx, idx, :len(sigma)] = sigma
        # Radiation
        for d_idx, dof in enumerate(dofs):
            res_rad = next((res for res in results_data if
                            res["type"] == "radiation" and
                            res["omega"] == omega and
                            res["dof"] == dof), None)
            if res_rad is not None:
                sigma = np.array(res_rad["sources"])
                SourcesR[w_idx, d_idx, :len(sigma)] = sigma*(-1j*omega)
                SourcesRb[w_idx, d_idx, :len(sigma)] = sigma
    return SourcesD, SourcesRb


def extract_potential_by_problem(results_data, unique_omegas, unique_betas, dofs, L, Ntheta, Nz):
    Nw = len(unique_omegas)
    Nb = len(unique_betas)
    Ndof = len(dofs)
    PotentialD = np.zeros((Nw, Nb*(L+1), Ntheta*Nz), dtype=complex)
    PotentialR = np.zeros((Nw, Ndof, Ntheta*Nz), dtype=complex)
    g=9.81
    for w_idx, omega in enumerate(unique_omegas):
        # Diffraction
        for l in range(L + 1):                # modes verticaux
            for b_idx, beta in enumerate(unique_betas):  # modes angulaires m
                res_diff = next((res for res in results_data if
                                    res["type"] == "diffraction" and
                                    res["omega"] == omega and
                                    res["beta"] == beta and
                                    res.get("modeL", 0) == l), None)
                if res_diff is not None:
                    phi = np.array(res_diff["cylsurface"])
                    # Calcul de l’indice dans la deuxième dimension : décalage selon l et m
                    idx = l * Nb + b_idx
                    PotentialD[w_idx, idx, :len(phi)] = phi*g/(1j*omega)
        # Radiation
        for d_idx, dof in enumerate(dofs):
            res_rad = next((res for res in results_data if
                            res["type"] == "radiation" and
                            res["omega"] == omega and
                            res["dof"] == dof), None)
            if res_rad is not None:
                phi = np.array(res_rad["cylsurface"])
                PotentialR[w_idx, d_idx, :len(phi)] = phi*g/(1j*omega)*(-1j*omega)
                # PotentialR[w_idx, d_idx, :len(phi)] = phi
    PotentialR=PotentialR.reshape(Nw, Ndof, Nz, Ntheta)
    PotentialD=PotentialD.reshape(Nw, Nb, Nz, Ntheta)
    return PotentialD, PotentialR



def Read_HydroCoefs(results_folder, unique_omegas, unique_betas, dofs):
    """
    - madd_array[w, i, j]
    - crad_array[w, i, j]
    - fex_array[w, b, i]
    """
    
    fe_abs_file = os.path.join(results_folder, "Capytaine_Fe_abs.dat")
    fe_phase_file = os.path.join(results_folder, "Capytaine_Fe_phase.dat")
    madd_file = os.path.join(results_folder, "Capytaine_Madd.dat")
    crad_file = os.path.join(results_folder, "Capytaine_Crad.dat")
	
    Nw = len(unique_omegas)
    Nb = len(unique_betas)
    Ndof = len(dofs)

    # --- FEX ---
    fex_array = np.zeros((Nw, Nb, Ndof), dtype=complex)
    with open(fe_abs_file, "r") as f_abs, open(fe_phase_file, "r") as f_phase:
        f_abs.readline(); f_phase.readline()  # Skip headers
        for line_abs, line_phase in zip(f_abs, f_phase):
            tokens_abs = line_abs.strip().split()
            tokens_phase = line_phase.strip().split()
            beta = float(tokens_abs[0])
            omega = float(tokens_abs[1])
            Fe_abs = list(map(float, tokens_abs[2:]))
            Fe_phase = list(map(float, tokens_phase[2:]))
            complex_vals = [a * np.exp(1j * p) for a, p in zip(Fe_abs, Fe_phase)]

            w_idx_array = np.where(np.isclose(unique_omegas, omega))[0]
            b_idx_array = np.where(np.isclose(unique_betas, beta))[0]

            if len(w_idx_array) == 0 or len(b_idx_array) == 0:
                print(f"[WARNING] Omega={omega:.6f} or beta={beta:.6f} not found in input lists.")
                continue

            w_idx = w_idx_array[0]
            b_idx = b_idx_array[0]
            fex_array[w_idx, b_idx, :] = complex_vals

    # --- MADD ---
    madd_array = np.zeros((Nw, Ndof, Ndof))
    with open(madd_file, "r") as f:
        f.readline()  # skip header

        while True:
            lines = [f.readline() for _ in range(Ndof)]
            if not lines[0]:
                break  # fin du fichier

            omega = float(lines[0].strip().split()[0])
            mat = []
            for line in lines:
                tokens = line.strip().split()
                coeffs = list(map(float, tokens[1:]))  # on saute omega
                mat.append(coeffs)

            mat = np.array(mat)
            if mat.shape != (Ndof, Ndof):
                raise ValueError(f"Madd matrix has unexpected shape: {mat.shape}")

            w_idx_array = np.where(np.isclose(unique_omegas, omega, rtol=1e-3))[0]
            if len(w_idx_array) == 0:
                print(unique_omegas)
                print(f"[WARNING] Omega={omega:.6f} not found in unique_omegas.")
                continue
            w_idx = w_idx_array[0]
            madd_array[w_idx, :, :] = mat


    # --- CRAD ---
    crad_array = np.zeros((Nw, Ndof, Ndof))
    with open(crad_file, "r") as f:
        f.readline()
        while True:
            lines = [f.readline() for _ in range(Ndof)]
            if not lines[0]:
                break  # fin du fichier

            omega = float(lines[0].strip().split()[0])
            mat = []
            for line in lines:
                tokens = line.strip().split()
                coeffs = list(map(float, tokens[1:]))  # on saute omega
                mat.append(coeffs)

            mat = np.array(mat)
            if mat.shape != (Ndof, Ndof):
                raise ValueError(f"Madd matrix has unexpected shape: {mat.shape}")

            w_idx_array = np.where(np.isclose(unique_omegas, omega, rtol=1e-3))[0]
            if len(w_idx_array) == 0:
                print(unique_omegas)
                print(f"[WARNING] Omega={omega:.6f} not found in unique_omegas.")
                continue
            w_idx = w_idx_array[0]
            crad_array[w_idx, :, :] = mat
            
    return madd_array, crad_array, fex_array


def Read_Potentials(unique_omegas, unique_betas, dofs, Ntheta, Nz, folder="./cylsurface"):
    """
    works for L = 0 ONLY
   
    """
    Nw = len(unique_omegas)
    Nb = len(unique_betas)
    Ndof = len(dofs)
    g = 9.81

    PotentialD = np.zeros((Nw, Nb, Ntheta * Nz), dtype=complex)
    PotentialR = np.zeros((Nw, Ndof, Ntheta * Nz), dtype=complex)

    pb_index = 0  

    for w_idx, omega in enumerate(unique_omegas):
        # Diffraction problems (one per beta)
        for b_idx, beta in enumerate(unique_betas):
            filename = os.path.join(folder, f"cylsurface.{pb_index+1:05d}.dat")
            pb_index += 1

            if not os.path.exists(filename):
                print(f"[WARNING] Fichier manquant : {filename}")
                continue

            try:
                data = np.loadtxt(filename)
                re, im = data[:, -2], data[:, -1]
                phi = re + 1j * im
                PotentialD[w_idx, b_idx, :len(phi)] = phi * g / (1j * omega)
            except Exception as e:
                print(f"[ERROR] Problème de lecture {filename} : {e}")

        # Radiation problems (one per dof)
        for d_idx, dof in enumerate(dofs):
            filename = os.path.join(folder, f"cylsurface.{pb_index+1:05d}.dat")
            pb_index += 1

            if not os.path.exists(filename):
                print(f"[WARNING] Fichier manquant : {filename}")
                continue

            try:
                data = np.loadtxt(filename)
                re, im = data[:, -2], data[:, -1]
                phi = re + 1j * im
                PotentialR[w_idx, d_idx, :len(phi)] = phi* g / (1j * omega)*(-1j*omega)
                # PotentialR[w_idx, d_idx, :len(phi)] = phi
            except Exception as e:
                print(f"[ERROR] Problème de lecture {filename} : {e}")

    PotentialD = PotentialD.reshape(Nw, Nb, Nz, Ntheta)
    PotentialR = PotentialR.reshape(Nw, Ndof, Nz, Ntheta)

    return PotentialD, PotentialR


def Read_Sources_From_Files(unique_omegas, unique_betas, dofs, L, Npanels, folder="./sources"):
    """
    Lit les fichiers sources.XXXXX.dat et reconstruit :
    - SourcesD : [Nw, (L+1)*Nb, Npanels]
    - SourcesRb : [Nw, Ndof, Npanels]
    - SourcesR : [-iω * SourcesRb]
    """
    Nw = len(unique_omegas)
    Nb = len(unique_betas)
    Ndof = len(dofs)

    SourcesD = np.zeros((Nw, (L + 1) * Nb, Npanels), dtype=complex)
    SourcesRb = np.zeros((Nw, Ndof, Npanels), dtype=complex)
    SourcesR = np.zeros((Nw, Ndof, Npanels), dtype=complex)

    pb_index = 0 

    for w_idx, omega in enumerate(unique_omegas):
        for l in range(L + 1):
            for b_idx, beta in enumerate(unique_betas):
                filename = os.path.join(folder, f"sources.{pb_index+1:05d}.dat")
                pb_index += 1

                if not os.path.exists(filename):
                    print(f"[WARNING] Fichier manquant : {filename}")
                    continue

                try:
                    data = np.loadtxt(filename)
                    re, im = data[:, 0], data[:, 1]
                    sigma = re + 1j * im
                    idx = l * Nb + b_idx  # index combiné (l, β)
                    SourcesD[w_idx, idx, :len(sigma)] = sigma
                except Exception as e:
                    print(f"[ERROR] Lecture {filename} : {e}")

        for d_idx, dof in enumerate(dofs):
            filename = os.path.join(folder, f"sources.{pb_index+1:05d}.dat")
            pb_index += 1

            if not os.path.exists(filename):
                print(f"[WARNING] Fichier manquant : {filename}")
                continue

            try:
                data = np.loadtxt(filename)
                re, im = data[:, 0], data[:, 1]
                sigma = re + 1j * im
                SourcesRb[w_idx, d_idx, :len(sigma)] = sigma
                SourcesR[w_idx, d_idx, :len(sigma)] = sigma * (-1j * omega)
            except Exception as e:
                print(f"[ERROR] Lecture {filename} : {e}")

    return SourcesD, SourcesRb

def wavenumber(folder):
    filename = os.path.join(folder, "Data_waves.dat")
    data=np.loadtxt(filename, skiprows=1)
    return data[:,1]

def WriteComputeTime(results_folder, time):
    os.makedirs(results_folder, exist_ok=True)

    # Préparation des fichiers
    file = os.path.join(results_folder, "Time.txt")
    with open(file, "w") as f:
        f.write(f"Computation Time (s) : {time:.4f}\n")



# Imports
import numpy as np
# from toolbox import mesh
from utils import configuration as CIT
from utils import output
from utils import plot
from utils import BEMoutput as out
import utils.files as UT
import time

# import utils.files as UT
import capytaine as cpt
from capytaine.post_pro import rao
from capytaine.bem.airy_waves import froude_krylov_force
from capytaine.bem.cylindrical_waves import froude_krylov_force_cyl
from capytaine.interaction_theory.BodyMesh import get_panel_centers_cylindrical, OCmesh
from capytaine.interaction_theory.Body import Body
from capytaine.interaction_theory.MultiBody import MultiBody as MB
import os
import json

project_data = json.load(open('capytaine_IT_projet.json', 'r'))

#################################################################
#################################################################
################### SINGLE BODY HYDRODYNAMICS ###################
#################################################################
#################################################################

## inputs to run nemoh
nemoh_def = project_data['nemoh_project']
dirwecs = nemoh_def['folder']
nemohdynamics = os.path.join(os.getcwd(), 'BEM', 'nemoh')
MeshFolder=os.path.join(dirwecs, 'MESH')
meshfile = os.path.join(MeshFolder, nemoh_def['mesh_filename'])
run_nemoh = nemoh_def['run_bem']
dire = os.path.join(dirwecs, nemoh_def['name'])
fdat = os.path.join(dire, f"{nemoh_def['mesh_filename'][:-4]}.dat")
convention = 'N'
results = os.path.join(dire, "resultsBEM")
outputBEM=True
resultsIT = os.path.join(dire, "resultsIT")
results_h5 = os.path.join(dire, "Farm_DMM.h5")
Motion = os.path.join(dire, "Motion")
PlotFolder = os.path.join(dire, "Plots")
mesh_ = os.path.join(dire, "mesh")
picklefile = os.path.join(dire, f"{nemoh_def['name']}_results.p") 
picklefile_E = os.path.join(dire, f"{nemoh_def['name']}_results_E.p") 


body_def = project_data['single_body_definition']
water_depth = body_def['water_depth'] # meters

min_value = body_def['spectral_param']['min']
max_value= body_def['spectral_param']['max']
format = body_def['spectral_param']['format']

if format == 'f':
    min_w = min_value*2*np.pi
    max_w = max_value*2*np.pi
elif format == 'T':
    min_w = 2*np.pi/min_value
    max_w = 2*np.pi/max_value
elif format == 'k':
    min_w = np.sqrt(9.81 * min_value * np.tanh(min_value * water_depth)) 
    max_w = np.sqrt(9.81 * max_value * np.tanh(max_value * water_depth)) 
elif format == 'L':
    min_w = np.sqrt(9.81 * (2*np.pi/min_value) * np.tanh((2*np.pi/min_value) * water_depth)) 
    max_w = np.sqrt(9.81 * (2*np.pi/max_value) * np.tanh((2*np.pi/max_value) * water_depth)) 
else : #'w'
    min_w = min_value
    max_w = max_value
if min_w>max_w :
        min_w, max_w=max_w, min_w
omega = np.linspace(min_w, max_w,
                    body_def['spectral_param']['number']) # rad/s 

if body_def['direction']['format'] == 'Nbeta':
    Nbeta=body_def['direction']['number']
elif body_def['direction']['format'] == 'M':
    Nbeta=2*body_def['direction']['number']+1
else :
    print("error : type not known, Nbeta or M")
directions = np.linspace(0.0, 360.0, Nbeta, endpoint=False) # degrees
directions *= np.pi/180.0 # radians                  

#################################################################
#################################################################
################### BEM with Capytaine ##############################
#################################################################
#################################################################
cpt.set_logging('WARNING')
geo = cpt.load_mesh(meshfile, file_format='nemoh')
# if nemoh_def["symmetry"]==1:
#     geo = cpt.ReflectionSymmetricMesh(geo, cpt.xOz_Plane, name=f"{nemoh_def['mesh_filename'][:-4]}")

if len(body_def['DOF']) == 6:
    # On crée directement avec rigid_body_dofs
    body = cpt.FloatingBody(
        mesh=geo,
        dofs=cpt.rigid_body_dofs(rotation_center=body_def['rotation_center']),
        center_of_mass=(0, 0, 0)
    )
    hydrostatics = body.compute_hydrostatics(rho=1000.0)

else:
    # On crée un corps vide et ajoute les DOFs une à une
    body = cpt.FloatingBody(mesh=geo, center_of_mass=(0, 0, 0))
    for dof in body_def['DOF']:
        if dof in ['Surge', 'Sway', 'Heave']:
            body.add_translation_dof(name=dof)
        elif dof in ['Roll', 'Pitch', 'Yaw']:
            body.add_rotation_dof(name=dof)
        else:
            raise ValueError(f"DOF inconnue : {dof}")

# print("Degrés de liberté du corps :", list(body.dofs.keys()))
vect_dof=list(body.dofs)
solver = cpt.BEMSolver()



print("------------------ BEM resolution with Capytaine -------------")
print("name : ", nemoh_def['name'])
print("Nw : ", len(omega))
print("Nbeta isolated : ", len(directions))
print("DOF : ", len(vect_dof), vect_dof)
print("Depth : ", water_depth)

UT.nemoh_structure(dire)
farm = project_data['farm_definition']
L=farm["Ne_modes"]
M=int((len(directions)-1)/2)
modeM=np.arange(-M, M+1)
modeL=np.arange(0, L+1)
if L>0 :
    print("-- with evanescent waves (Ne modes =", L, ")")
if body_def['sources']==1:
    Incident_potential="cylindrical"
else :
    Incident_potential="planar"
    OC=body_def['cylinder']
    OuterCylinder=OCmesh(OC['radius'], OC['azimuthal_discretization'], OC['depth_discretization'], water_depth)
results_data = []
wavenumber=np.zeros(len(omega))
Np=len(omega)*(len(directions)*(L+1)+len(vect_dof))

if run_nemoh : 
    start_time = time.time()
    num = 1
    for i, w in enumerate(omega):
        print("problem n°", num, "/",Np)

        for e in modeL :
            for k, angle in enumerate(directions) :
            # --- DIFFRACTION ---
                pbD = cpt.DiffractionProblem(body=body, wave_direction=angle, 
                            omega=w, water_depth=water_depth, BC = Incident_potential, modeM=modeM[k], modeL=e)
                resultD = solver.solve(pbD)
                if body_def['sources']==1:
                    fk_force = froude_krylov_force_cyl(pbD)
                else :
                    fk_force = froude_krylov_force(pbD)

                excitation_force = {dof: fk_force.get(dof, 0.0) + resultD.forces.get(dof, 0.0) for dof in resultD.forces}
                # excitation_force = resultD.forces  # dict par DOF
                sources_diff = resultD.sources    
                if body_def['sources']==1:
                    cylsurface = None
                else :
                    cylsurface = solver.compute_cylsurface_elevation(OuterCylinder, resultD)

                results_data.append({
                "type": "diffraction",
                "omega": w,
                "beta": angle,
                "dof": None,
                "modeM" : modeM[k],
                "modeL" : e,
                "excitation_force": excitation_force,
                "sources": sources_diff,
                "cylsurface": cylsurface,
                "added_mass": None,
                "damping": None
                })
                num += 1
        # --- RADIATION ---
        for dof in vect_dof:
            pbR = cpt.RadiationProblem(body=body, radiating_dof=dof, omega=w, water_depth=water_depth)
            resultR = solver.solve(pbR)

            added_mass = resultR.added_mass
            damping = resultR.radiation_damping
            sources_rad = resultR.sources
            if body_def['sources']==1:
                cylsurfaceR = None
            else :
                cylsurfaceR = solver.compute_cylsurface_elevation(OuterCylinder, resultR)

            results_data.append({
            "type": "radiation",
            "omega": w,
            "beta": None,
            "dof": dof,
            "modeM" : None,
            "modeL" : None,
            "excitation_force": None,
            "sources": sources_rad,
            "cylsurface": cylsurfaceR,
            "added_mass": added_mass,
            "damping": damping
            })

            num += 1
        wavenumber[i]=pbR.wavenumber
    end_time = time.time()
    out.WriteComputeTime(results, end_time-start_time)
    out.write_params_file(results_data, water_depth, output_file=os.path.join(dire, "params.dat"))    
    if not os.path.exists(results):
        os.makedirs(results)
    if body_def['sources']==0:
        out.WriteCylsurfaceByProblem(results_data, OuterCylinder, results_folder=results)
    else:
        out.WriteSourcesByProblem(results_data, body, results_folder=results)
    out.WriteCapytaineDataFromResults(results_data, results_folder=results, dofs=vect_dof)

    Madd_iso, Crad_iso, Fex_iso = out.extract_hydrodynamic_quantities_with_beta(results_data, vect_dof, omega, directions, L)
    if body_def['sources']==1:
        centers, areas=get_panel_centers_cylindrical(meshfile)
        SourcesS, SourcesR = out.extract_sources_by_problem(results_data, omega, directions, vect_dof, L, len(areas))
    else : 
        phiS, phiR = out.extract_potential_by_problem(results_data, omega, directions, vect_dof, L, OC['azimuthal_discretization'], OC['depth_discretization'])
else :
    # print("to do")
    # out.CheckData() here or in Body 
    wavenumber = out.wavenumber(results)
    Madd_iso, Crad_iso, Fex_iso = out.Read_HydroCoefs(results, omega, directions, vect_dof)
    if body_def['sources']==1:
        centers, areas=get_panel_centers_cylindrical(meshfile)
        SourcesS, SourcesR = out.Read_Sources_From_Files(omega, directions, vect_dof, L, len(areas), results)
    else : 
        phiS, phiR = out.Read_Potentials(omega, directions, vect_dof, OC['azimuthal_discretization'], OC['depth_discretization'], results)
    
if body_def['sources']==1:
    WEC = Body(omega, directions, water_depth, L, wavenumber, Fex_iso, True)
    WEC.Transfers(SourcesS, SourcesR, (centers, areas))
else : 
    WEC = Body(omega, directions, water_depth, L, wavenumber, Fex_iso, False)
    # OC_th = np.array([OuterCylinder[j][1] for j in range(OC['azimuthal_discretization'])])  
    OC_th = np.linspace(0, 2*np.pi, OC['azimuthal_discretization'], endpoint = False)   
    OC_z = np.array([OuterCylinder[i * OC['azimuthal_discretization']][2] for i in range(OC['depth_discretization'])])
    WEC.Transfers(phiS, phiR, (OC['radius'], OC_th,OC_z))

WEC.Write(results, L, body_def['sources'])
#################################################################
#################################################################
################### DIRECT MATRIX METHOD ########################
#################################################################
#################################################################
print("------------------ Interaction Theory -------------")

## inputs to generate array
betas = farm['wave_directions']
directionMB = np.linspace(betas['min'],
                        betas['max'],
                        betas['number']) # must be converted to degrees WHY ?? NO
print("beta=", directionMB)
if betas['format'] == 'DEG':
    # directionMB *= 180.0/np.pi  #RAD to DEG
    directionMB *= np.pi/180.0   #DEG to RAD
print("Nbeta system : ", len(directionMB))


param_distance=2
limite=600
# radius_barge=6.36
diameter= 2*body_def['cylinder']['radius']

while param_distance<=limite : 
    if farm['Configuration']=="logd":
        N_bodies=farm["N_bodies"]
        print("Nb : ", N_bodies)
        if N_bodies>1 :
            distance = param_distance
            print(f"Distance between the two bodies center : {distance+diameter}")
            print(f"Distance between the two bodies : {distance}")
        else :
            distance=0
        print("Type configuration", farm["Type"])
        coord = CIT.CreateConfig(N_bodies, diameter+distance, farm["Type"])
        limite=600
    elif farm['Configuration']=="d/a":
        N_bodies=farm["N_bodies"]
        print("Nb : ", N_bodies)
        radius_a=farm["body_dim"]
        if N_bodies>1 :
            distance = param_distance
            print(f"Distance between the two bodies center : {radius_a*distance}")
            print(f"Distance / a : {distance}")
        else :
            distance=0
        print("Type configuration", farm["Type"])
        coord = CIT.CreateConfig(N_bodies, radius_a*distance, farm["Type"])
        limite=10
    else :
        layout = farm['layout']
        print("Nb : ", len(layout))
        N_bodies=len(layout)
        coord = np.zeros((len(layout), 2))
        labels = []
        for ib, bodyP in enumerate(layout):
            coord[ib, :] = bodyP["position"]
            labels.append(bodyP["name"])
        if len(coord)>1 :
            distance = np.linalg.norm(coord[1] - coord[0])
            print(f"Distance between the two bodies center : {distance}")
            # distance = distance - farm["body_dim"]*2
            a=farm["body_dim"]
            print(f"Distance between the two bodies : {distance- a*2}")
        else :
            distance=0
        param_distance=limite+1

    WECArr = MB(directionMB, WEC, len(vect_dof), Madd_iso, Crad_iso)
    WECArr.Scattering(coord, L) #  WECArr.Fex.shape (Num freq, Num dir, Num dof)
    WECArr.Radiation(coord, L) # Available: WECArr.Madd and WECArr.Crad, which stand for added mass and radiation damping, respectively: WECArr.Madd.shape (Num freq, Num dof, Num dof)
    
    if farm['RAO'] :
        print("-- RAO calculation activated")
        inertia_matrix = body.build_inertia_matrix_from_translated_bodies(coord) 
        # print("M", inertia_matrix)
        hydrostatic_stiffness = body.build_hydrostatic_stiffness_from_translated_bodies(coord)
        # print("K", hydrostatic_stiffness)
        if not os.path.exists(Motion):
            os.makedirs(Motion)
        WECArr.RAO(os.path.join(dire, "Mechanics"), inertia_matrix, hydrostatic_stiffness) 

    if farm['Free_surface']['Nx'] > 0 :
        print("-- Free surface calculation activated")
        if not os.path.exists(Motion):
            os.makedirs(Motion)
        WECArr.FreeSurface(coord, farm['Free_surface']['Nx'], farm['Free_surface']['Ny'], 
                           farm['Free_surface']['Lx'], farm['Free_surface']['Ly'], farm["Ne_modes"]) 

#################################################################
#################################################################
######################### OUTPUT ################################
#################################################################
#################################################################
    if not os.path.exists(resultsIT):
        os.makedirs(resultsIT)
    out=output.WriteData(WECArr, directionMB, farm, distance, resultsIT, results_h5, Motion, body_def['sources'])
    if len(farm['Plot_DOF'])>0 :
        if not os.path.exists(PlotFolder):
            os.makedirs(PlotFolder)
        out=plot.PlotData(WECArr, farm, distance, directionMB, farm['Plot_DOF'], PlotFolder)
    
    if farm['Configuration']=="d/a":
        param_distance=param_distance+1
    else :
        param_distance=param_distance*2


assert len(WECArr.Fex[:, 0, 0]) == len(WECArr.Madd[:, 0, 0]) == len(WECArr.period)
assert len(WECArr.Fex[0, 0, :]) == len(WECArr.Madd[0, :, 0]) == len(WECArr.Madd[0, 0, :]) == N_bodies*len(vect_dof)
assert len(WECArr.Fex[0, :, 0]) == len(directionMB)

## NOTE: Notice you only need the instance WEC to run the DIRECT MATRIX METHOD.
## Once you've got WEC you can skip the part SINGLE BODY HYDRODYNAMICS and go straight to DIRECT MATRIX METHOD for different array configurations (coord) or different wave headings (directionMB).
# plt.show()
print('end of script')

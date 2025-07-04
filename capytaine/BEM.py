

# Imports
import numpy as np
# from toolbox import mesh
from utils import configuration as CIT
from utils import output
from utils import plot
from utils import BEMoutput as out
import utils.files as UT

# import utils.files as UT
import capytaine as cpt
from capytaine.post_pro import rao
from capytaine.bem.airy_waves import froude_krylov_force
from capytaine.bem.cylindrical_waves import froude_krylov_force_cyl
import os
import json

project_data = json.load(open('capytaine_BEM_projet.json', 'r'))

#################################################################
#################################################################
################### SINGLE BODY HYDRODYNAMICS ###################
#################################################################
#################################################################

## inputs to run nemoh
nemoh_def = project_data['nemoh_project']
dirwecs = nemoh_def['folder']
MeshFolder=os.path.join(dirwecs, 'MESH')
meshfile = os.path.join(MeshFolder, nemoh_def['mesh_filename'])
dire = os.path.join(dirwecs, nemoh_def['name'])
fdat = os.path.join(dire, f"{nemoh_def['mesh_filename'][:-4]}.dat")
results = os.path.join(dire, "resultsBEM")
Motion = os.path.join(dire, "Motion")
PlotFolder = os.path.join(dire, "Plots")
mesh_ = os.path.join(dire, "mesh")


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

directions = np.linspace(body_def['direction']['min'],
                            body_def['direction']['max'],
                            body_def['direction']['number'],
                            endpoint=False) # degrees
if body_def['direction']['format'] == 'DEG':
    directions *= np.pi/180.0


#################################################################
#################################################################
################### BEM with Capytaine ##############################
#################################################################
#################################################################

print("------------------ BEM resolution with Capytaine -------------")
print("name : ", nemoh_def['name'])
print("Nw : ", len(omega))
print("Nbeta isolated : ", len(directions))
cpt.set_logging('WARNING')
barge = cpt.load_mesh(meshfile, file_format='nemoh')
solver = cpt.BEMSolver()
if body_def['sources']==1:
    Incident_potential="cylindrical"
    S="S_"
else :
    Incident_potential="planar"
    S=''
param_distance=1
limite=1023

diameter= 2*body_def['cylinder']['radius']
farm = project_data['farm_definition']
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
        results_d = os.path.join(results, S+farm["Type"]+f"_d{distance}")

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
        for ib, body in enumerate(layout):
            coord[ib, :] = body["position"]
            labels.append(body["name"])
        if len(coord)>1 :
            distance = np.linalg.norm(coord[1] - coord[0])
            print(f"Distance between the two bodies center : {distance}")
            # distance = distance - diameter
            print(f"Distance between the two bodies : {distance}")
        else :
            distance=0
        results_d = os.path.join(results, S+f"Nb{N_bodies}_d{distance}")
        param_distance=limite+1
 
    if farm['Configuration']=="d/a":
        param_distance=param_distance+1
    else :
        param_distance=param_distance*2
    body = cpt.FloatingBody(mesh=barge,
                        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
                        center_of_mass=(0, 0, 0))
    # locations = np.array([[0.0, 0.0], [10.0, 0.0]])
    all_bodies = body.assemble_arbitrary_array(coord)
    vect_dof=list(all_bodies.dofs)
    print("DOF : ", len(vect_dof))
    UT.nemoh_structure(dire)
    M=int((len(directions)-1)/2)
    modeM=np.arange(-M, M+1)
    
    results_data = []
    wavenumber=np.zeros(len(omega))
    Np=len(omega)*(len(directions)+len(vect_dof))
    num = 1
    for i, w in enumerate(omega):
        print("problem n°", num, "/",Np)

        for k, angle in enumerate(directions) :
        # --- DIFFRACTION ---
            pbD = cpt.DiffractionProblem(body=all_bodies, wave_direction=angle, 
                            omega=w, water_depth=water_depth, BC=Incident_potential, modeM=modeM[k])
            resultD = solver.solve(pbD)
            
            if body_def['sources']==1:
                fk_force = froude_krylov_force_cyl(pbD)
            else :
                fk_force = froude_krylov_force(pbD)
            excitation_force = {
            dof: fk_force.get(dof, 0.0) + resultD.forces.get(dof, 0.0)
            for dof in resultD.forces
            }
            # excitation_force = resultD.forces  # dict par DOF
            sources_diff = resultD.sources    

            results_data.append({
                "type": "diffraction",
                "omega": w,
                "beta": angle,
                "dof": None,
                "modeM" : modeM[k],
                "modeL" : None,
                "excitation_force": excitation_force,
                "sources": sources_diff,
                "added_mass": None,
                "damping": None
            })
            num += 1
        # --- RADIATION ---
        for dof in vect_dof:
            pbR = cpt.RadiationProblem(body=all_bodies, radiating_dof=dof, omega=w, water_depth=water_depth)
            resultR = solver.solve(pbR)

            added_mass = resultR.added_mass
            damping = resultR.radiation_damping
            sources_rad = resultR.sources

            results_data.append({
                "type": "radiation",
                "omega": w,
                "beta": None,
                "dof": dof,
                "modeM" : None,
                "modeL" : None,
                "excitation_force": None,
                "sources": sources_rad,
                "added_mass": added_mass,
                "damping": damping
            })

            num += 1
        wavenumber[i]=pbR.wavenumber
    out.write_params_file(results_data, output_file=os.path.join(dire, "params.dat"))    
    if not os.path.exists(results):
        os.makedirs(results)
    if not os.path.exists(results_d):
        os.makedirs(results_d)
    print(results_d)
    out.WriteCapytaineDataFromResults(results_data, results_folder=results_d, dofs=vect_dof)
    out.WriteSourcesByProblem(results_data, all_bodies, results_folder=results_d)

print('end of script')

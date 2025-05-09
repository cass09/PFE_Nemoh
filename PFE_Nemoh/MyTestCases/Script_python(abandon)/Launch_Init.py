# Imports
import numpy as np
from tools_python import mesh
import tools_python.files as UT
import json
import matplotlib.pyplot as plt
import h5py

import os
import json

project_data = json.load(open('Nemoh_project.json', 'r'))

#################################################################
#################################################################
################### SINGLE BODY HYDRODYNAMICS ###################
#################################################################
#################################################################

## inputs to run nemoh
nemoh_def = project_data['nemoh_project']
dirwecs = nemoh_def['folder']
# nemohdynamics = os.path.join(os.getcwd(), 'BEM', 'nemoh')
dirmesh=os.path.join(dirwecs, "MESH")
fdat = os.path.join(dirmesh, nemoh_def['mesh_filename'])
run_nemoh = nemoh_def['run_bem']
dire = os.path.join(dirwecs, nemoh_def['name'])
new_fdat = os.path.join(dire, nemoh_def['mesh_filename'])

results = os.path.join(dire, "results")
results_h5 = os.path.join(dire, "Farm_DMM.h5")

mesh_ = os.path.join(dire, "mesh")
picklefile = os.path.join(dire, f"{nemoh_def['name']}_results.p") 


body_def = project_data['single_body_definition']
freqs = np.linspace(body_def['frequency']['min'],
                    body_def['frequency']['max'],
                    body_def['frequency']['number']) # rad/s
    
if body_def['frequency']['format'] == 'HZ':
    freqs *= 2*np.pi

directions = np.linspace(body_def['direction']['min'],
                            body_def['direction']['max'],
                            body_def['direction']['number'],
                            endpoint=True) # degrees

if body_def['direction']['format'] == 'RAD':
    directions *= 180.0/np.pi

Env = project_data['environment']
Depth=Env['water_depth']
Rho=Env["rho"]
G=Env["Gravity"]
REF=Env["REF"]


# print("name : ", nemoh_def['name'])
# print("Nw : ", len(freqs))
# print("Nbeta isolated : ", len(directions))
# print("Nmodes : ", len(body_def['modes']))


if run_nemoh:    
    UT.nemoh_structure(dire)
    BEM=1
else :
    BEM=0
    
translation = np.zeros(3)
rotation = np.zeros(3)
geom = mesh.readDAT(fdat, translation, rotation)
if np.any(translation != np.zeros(3)) or np.any(rotation != np.zeros(3)):
    modify_meshfile=True
else :
    modify_meshfile=False
geom.save_mesh(fdat, new_fdat, modify_meshfile)

## info
mode = body_def['modes']
modes = ([mode, mode],)

Cylsurface = (body_def['cylinder']['radius'],
                body_def['cylinder']['azimuthal_discretization'],
                body_def['cylinder']['depth_discretization']) # (radius of the cylinder enclosing the body, azimuthal discretization, depth discretization)

meshes = ([geom.dirDAT, geom.Nnodes, geom.Npanels], )

## generate inputs for running nemoh
freqs_nemoh = (len(freqs), freqs.min(), freqs.max())
directions_nemoh = (len(directions), directions.min(), directions.max())
# directory = UT.InputDynamics(dire, meshes, modes, freqs_nemoh, directions_nemoh, EnvData, Cylsurface)
    
    


## inputs to generate array
farm = project_data['farm_definition']
layout = farm['layout']
print("Nb : ", len(layout))
coord = np.zeros((len(layout), 2))
labels = []
for ib, body in enumerate(layout):
    coord[ib, :] = body["position"]
    labels.append(body["name"])

if len(coord)>1 :
    distance = np.linalg.norm(coord[1] - coord[0])
    print(f"Distance between the two bodies center : {distance}")
    distance = distance - 2*body_def['cylinder']['radius']
    print(f"Distance between the two bodies : {distance}")
else :
    distance=0


betas = farm['wave_directions']

directionMB = np.linspace(betas['min'],
                          betas['max'],
                          betas['number']) # must be converted to degrees
if betas['format'] == 'RAD':
    directionMB *= 180.0/np.pi
IT_directions = (len(directionMB), directionMB.min(), directionMB.max())
# print("Nbeta system : ", len(directionMB))


Post_Process=project_data["post_process"]
if Post_Process["IRF"]==False :
    IRF=0
else :
    IRF=1
dt=Post_Process["time_step"]
Tf=Post_Process["duration"]
if Post_Process["Pressure"]==False :
    P=0
else :
    P=1
Kochin=Post_Process["Kochin"]
if Kochin[0]["Active"]==False :
    K=[0, 0, 0, 0]
else :
    K=Kochin["params"]
FreeSurface=Post_Process["Free surface"]
if FreeSurface[0]["Active"]==False :
    FS=[0, 0, 0, 0]
else :
    FS=FreeSurface["params"]
if Post_Process["RAO"]==False :
    RAO=0
else :
    RAO=1
Uout=Post_Process["output_freq"]

Input_solver=project_data["solver"]
Gauss=Input_solver["N_Gauss"]
epsilon=Input_solver["eps_zmin"]
type_solver=Input_solver["solver_type"]
restart=Input_solver["restart"]
tol=Input_solver["tol"]
MaxIter=Input_solver["MaxIter"]



    # Shortcuts
_md = os.mkdir
_j = os.path.join


Nb = len(modes)

# make a mesh directory and a results directory
_md(_j(dire, 'mesh'))
_md(_j(dire, 'results'))

# Generate ID.dat
with open(_j(dire,'ID.dat'), 'w') as f :
    f.write('1\n')
    f.write('.')

# Generate input.txt
with open(_j(dire,'input_solver.txt'), 'w') as f :
    f.write('--- Calculation parameters --------------------------\n')
    f.write('0				! Solver selection')

# Generate Nemoh.cal
with open(_j(dire,'Nemoh.cal'), 'w') as f :
    f.write('--- Environment ------------------------------------------------------------------------------------------------------------------\n')
    f.write('{:.2f}				! RHO 			KG/M**3\n'.format(Rho))
    f.write('{:.2f}				! G			M/S**2\n'.format(G))
    f.write('{:.2f}				! WATER DEPTH			M (0 for infinite water depth)\n'.format(abs(Depth)))
    f.write('{:} {:}			    ! COORDINATES OF WAVE MEASUREMENT POINT         M	\n'.format(*REF))
    f.write('--- Description of floating bodies -----------------------------------------------------------------------------------------------\n')
    f.write('{:}				! Number of bodies\n'.format(Nb))
    for nb in range(Nb):
        f.write('--- Body {:} -------------------------------------------------------------------------------------------------------------------\n'.format(nb+1))
        f.write('{:}                ! Name of mesh file\n'.format(nemoh_def['mesh_filename']))
        f.write('{:} {:}			! Number of nodes and number of panels \n'.format(meshes[nb][1],meshes[nb][2]))
        Nd = len(modes[nb][0])
        f.write('{:}				! Number of generalized degrees of freedom\n'.format(Nd))
        for dof in modes[nb][0]:
            f.write('{:} {:} {:} {:} {:} {:} {:}		! DoF: 1 for translation (2 for roation), followed by direction vector (axis and center of rotation for rotation)\n'.format(*dof))
        Nf = len(modes[nb][1])
        f.write('{:}				! Number of generalized degrees of freedom\n'.format(Nf))
        for force in modes[nb][1]:
            f.write('{:} {:} {:} {:} {:} {:} {:}		! DoF: 1 for translation (2 for roation), followed by direction vector (axis and center of rotation for rotation)\n'.format(*force))
        f.write('0                                                    ! Number of additional lines\n')
    f.write('--- Load cases to be solved -------------------------------------------------------------------------------------------------------\n')
    f.write('{:} {:} {:}		! Number of wave frequencies, min, max (rad/s)\n'.format(*freqs_nemoh))
    f.write('{:} {:} {:}		! Number of wave direction, min, max (degrees)\n'.format(*directions_nemoh))
    f.write('--- Post processing ---------------------------------------------------------------------------------------------------------------\n')
    f.write('{:} {:} {:}		! IRF calculation (1 for yes), time step and duration\n'.format(IRF, dt, Tf))
    f.write('{:}				! Save pressure on body surface (1 for yes)\n'.format(P))
    f.write('{:} {:} {:}		! Kochin function calculation: nb of direction (0 for no), min, max (degrees)\n'.format(*K))
    f.write('{:} {:} {:} {:}	! Free surface vizualization: nb of points in x (0 for no), y and dimensions of domain in x and y directions\n'.format(*FS))
    f.write('{:}	            ! Response Amplitude Operator (RAO), 0 no calculation, 1 calculated -> Inertia.cal\n'.format(RAO))
    f.write('{:}	            ! output freq type, 1,2,3=[rad/s,Hz,s]\n'.format(Uout))
    f.write('--- Interaction Theory ---------------------------------------------------------------------------------------------------------------\n')
    f.write('1	                ! run Interaction Theory \n')
    f.write('{:}				! run BEM (1 for yes)\n'.format(BEM))        
    f.write('{:} {:} {:}		! Field points calculation: radius of the cylinder (0 for no calculations), number of points in theta and number of points in z\n'.format(*Cylsurface))
    f.write('{:}				! Number of bodies\n'.format(len(layout)))
    for nb in range(len(layout)):
        f.write('{:} {:}			! Coordinates body {:} \n'.format(*coord[nb], nb+1))
    f.write('{:} {:} {:}		! Number of wave direction, min, max (degrees)\n'.format(*IT_directions))
    f.write('--- QTF ---------------------------------------------------------------------------------------------------------------\n')
    f.write('0  				! QTF flag, 1 is calculated \n')
    f.write('\n')


# PB :  calcul Nb noeuds !
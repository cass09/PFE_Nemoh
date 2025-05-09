

# Imports
import numpy as np
from BEM import Nemoh
from semi_analytical.Body import Body
from semi_analytical.MultiBody import MultiBody as MB
from toolbox import mesh
from utils import configuration as CIT
from utils import output 
import utils.files as UT
import json
import matplotlib.pyplot as plt
import h5py
import cmath

import os
import json

project_data = json.load(open('project_definition.json', 'r'))

#################################################################
#################################################################
################### SINGLE BODY HYDRODYNAMICS ###################
#################################################################
#################################################################

## inputs to run nemoh
nemoh_def = project_data['nemoh_project']
dirwecs = nemoh_def['folder']
nemohdynamics = os.path.join(os.getcwd(), 'BEM', 'nemoh')
fgdf = os.path.join(dirwecs, nemoh_def['mesh_filename'])
run_nemoh = nemoh_def['run_bem']
dire = os.path.join(dirwecs, nemoh_def['name'])
fdat = os.path.join(dire, f"{nemoh_def['mesh_filename'][:-4]}.dat")
convention = 'N'
results = os.path.join(dire, "results")
results_h5 = os.path.join(dire, "Farm_DMM.h5")
Motion = os.path.join(dire, "Motion")
mesh_ = os.path.join(dire, "mesh")
picklefile = os.path.join(dire, f"{nemoh_def['name']}_results.p") 
picklefile_E = os.path.join(dire, f"{nemoh_def['name']}_results_E.p") 


body_def = project_data['single_body_definition']
freqs = np.linspace(body_def['frequency']['min'],
                    body_def['frequency']['max'],
                    body_def['frequency']['number']) # rad/s 
if body_def['frequency']['format'] == 'HZ':
    freqs *= 2*np.pi
directions = np.linspace(body_def['direction']['min'],
                            body_def['direction']['max'],
                            body_def['direction']['number'],
                            endpoint=False) # degrees
if body_def['direction']['format'] == 'RAD':
    directions *= 180.0/np.pi
depth = body_def['water_depth'] # meters

#################################################################
#################################################################
################### BEM with Nemoh ##############################
#################################################################
#################################################################

if run_nemoh:
    ## read gdf mesh file and create dat mesh file
    ## create the Nemoh folder structure or ask user for inputs
    print("------------------ Running Nemoh -----------------")
    UT.nemoh_structure(dire)
    
    if nemoh_def['mesh_type'] == "WAMIT" : 
        geom = mesh.readGDF(fgdf)
        geom.dat(fdat + '.dat')
    else : 
        geom = mesh.readDAT(fgdf)
        geom.dat(fdat + '.dat')
    
    mode = body_def['modes']
    modes = ([mode, mode],)
    
    FieldPoints = (body_def['cylinder']['radius'],
                   body_def['cylinder']['azimuthal_discretization'],
                   body_def['cylinder']['depth_discretization']) # (radius of the cylinder enclosing the body, azimuthal discretization, depth discretization)
    meshes = ([f"{nemoh_def['mesh_filename'][:-4]}.dat", geom.Nnodes, geom.Npanels], )
    ## generate inputs for running nemoh
    freqs_nemoh = (len(freqs), freqs.min(), freqs.max())
    directions_nemoh = (len(directions), directions.min(), directions.max())

    directory = Nemoh.InputDynamics(dire, meshes, modes, freqs_nemoh, 
                                    directions_nemoh, depth, FieldPoints)
    ## run nemoh
    status = Nemoh.RunDynamics(directory, nemohdynamics)


print("------------------ Interaction Theory -------------")
print("name : ", nemoh_def['name'])
print("Nw : ", len(freqs))
print("Nbeta isolated : ", len(directions))
print("DOF : ", len(body_def['modes']))

farm = project_data['farm_definition']
if farm['Evanescent'] :
    print("-- with evanescent waves (Ne modes =", farm["Ne_modes"], ")")
if farm['RAO'] :
    print("-- RAO calculation activated")

## generate wec (diffraction and force transfer matrices and radiation coefficents)
WEC = Body(freqs, directions*np.pi/180., depth, farm["Ne_modes"], convention)

if farm["Evanescent"] : 
    WEC.Transfers(dire, results, results, farm["Evanescent"], BEM='N', Tol= 1e-6)
    WEC.PickTransfers(picklefile_E, save=True)
    loaded_WEC = WEC.PickTransfers(picklefile_E, save=False)
else : 
    if not UT.file_exist(picklefile):
        WEC.Transfers(dire, results, results, farm["Evanescent"], BEM='N', Tol= 1e-6)
        WEC.PickTransfers(picklefile, save=True)
    
    loaded_WEC = WEC.PickTransfers(picklefile, save=False)
WEC.Write(results, farm["Evanescent"])


#################################################################
#################################################################
################### DIRECT MATRIX METHOD ########################
#################################################################
#################################################################

## inputs to generate array
betas = farm['wave_directions']
directionMB = np.linspace(betas['min'],
                        betas['max'],
                        betas['number']) # must be converted to degrees
if betas['format'] == 'RAD':
    directionMB *= 180.0/np.pi
print("Nbeta system : ", len(directionMB))


param_distance=1
limite=1023
# radius_barge=6.36
diameter= 2*body_def['cylinder']['radius']

while param_distance<=limite : 
    if farm['Configuration'] :
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
            distance = distance - diameter
            print(f"Distance between the two bodies : {distance}")
        else :
            distance=0
        param_distance=limite+1

    WECArr = MB(directionMB, WEC, cylamplitude=True)
    WECArr.Scattering(coord, farm["Evanescent"]) # Available: WECArr.Fex, which stands for excitation force: WECArr.Fex.shape (Num freq, Num dir, Num dof)
    WECArr.Radiation(coord, farm["Evanescent"]) # Available: WECArr.Madd and WECArr.Crad, which stand for added mass and radiation damping, respectively: WECArr.Madd.shape (Num freq, Num dof, Num dof)
    
    if farm['RAO'] :
        if not os.path.exists(Motion):
            os.makedirs(Motion)
        Mechanics = os.path.join(dire, "Mechanics")
        directoryH=Nemoh.InputHydro(dire, coord)
        statusH = Nemoh.RunHydro(directoryH, nemohdynamics)
        WECArr.RAO(Mechanics) 


#################################################################
#################################################################
######################### OUTPUT ################################
#################################################################
#################################################################
    out=output.WriteData(WECArr, directionMB, farm, distance, results, results_h5, Motion)
    param_distance=param_distance*2


assert len(WECArr.Fex[:, 0, 0]) == len(WECArr.Madd[:, 0, 0]) == len(WECArr.period)
assert len(WECArr.Fex[0, 0, :]) == len(WECArr.Madd[0, :, 0]) == len(WECArr.Madd[0, 0, :]) == N_bodies*len(body_def['modes'])
assert len(WECArr.Fex[0, :, 0]) == len(directionMB)


# fig, axs = plt.subplots(len(coord), 1)
# fig.suptitle('Excitation force coefficients')

# for ix, ax in enumerate(axs):
#     ax.plot(WECArr.period, np.abs(WECArr.Fex[:,0,ix]))
#     ax.set_title(labels[ix])
#     ax.set(xlabel='frequency, [rad/s]', ylabel='abs(Fex)')

# fig, axs = plt.subplots(len(coord), len(coord))
# fig.suptitle('Radiation force coefficients')

# for ii in range(len(coord)):
#     for jj in range(len(coord)):
#         axs[ii,jj].plot(WECArr.period, np.abs(WECArr.Madd[:,ii,jj]))
#         axs[ii,jj].plot(WECArr.period, np.abs(WECArr.Crad[:,ii,jj]))
#         ax.set_title(labels[ii])
#         ax.set(xlabel='frequency, [rad/s]', ylabel='Rad Damping')

## NOTE: Notice you only need the instance WEC to run the DIRECT MATRIX METHOD.
## Once you've got WEC you can skip the part SINGLE BODY HYDRODYNAMICS and go straight to DIRECT MATRIX METHOD for different array configurations (coord) or different wave headings (directionMB).
# plt.show()
print('end of script')
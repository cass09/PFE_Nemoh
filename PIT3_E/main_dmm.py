

# Imports
import numpy as np
from BEM import Nemoh
from semi_analytical.Body import Body
from semi_analytical.MultiBody import MultiBody as MB
from toolbox import mesh
from utils import configuration as CIT
from utils import output
from utils import plot
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
MeshFolder=os.path.join(dirwecs, 'MESH')
fgdf = os.path.join(MeshFolder, nemoh_def['mesh_filename'])
run_nemoh = nemoh_def['run_bem']
dire = os.path.join(dirwecs, nemoh_def['name'])
fdat = os.path.join(dire, f"{nemoh_def['mesh_filename'][:-4]}.dat")
convention = 'N'
results = os.path.join(dire, "results")
resultsIT = os.path.join(dire, "resultsIT")
results_h5 = os.path.join(dire, "Farm_DMM.h5")
Motion = os.path.join(dire, "Motion")
PlotFolder = os.path.join(dire, "Plots")
mesh_ = os.path.join(dire, "mesh")
picklefile = os.path.join(dire, f"{nemoh_def['name']}_results.p") 
picklefile_E = os.path.join(dire, f"{nemoh_def['name']}_results_E.p") 


body_def = project_data['single_body_definition']
depth = body_def['water_depth'] # meters

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
    min_w = np.sqrt(9.81 * min_value * np.tanh(min_value * depth)) 
    max_w = np.sqrt(9.81 * max_value * np.tanh(max_value * depth)) 
elif format == 'L':
    min_w = np.sqrt(9.81 * (2*np.pi/min_value) * np.tanh((2*np.pi/min_value) * depth)) 
    max_w = np.sqrt(9.81 * (2*np.pi/max_value) * np.tanh((2*np.pi/max_value) * depth)) 
else : #'w'
    min_w = min_value
    max_w = max_value
if min_w>max_w :
        min_w, max_w=max_w, min_w
freqs = np.linspace(min_w, max_w,
                    body_def['spectral_param']['number']) # rad/s 

directions = np.linspace(body_def['direction']['min'],
                            body_def['direction']['max'],
                            body_def['direction']['number'],
                            endpoint=False) # degrees
if body_def['direction']['format'] == 'RAD':
    directions *= 180.0/np.pi


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
                                    directions_nemoh, depth, FieldPoints, body_def['solver'])
    ## run nemoh
    status = Nemoh.RunDynamics(directory, nemohdynamics)


print("------------------ Interaction Theory -------------")
print("name : ", nemoh_def['name'])
print("Nw : ", len(freqs))
print("Nbeta isolated : ", len(directions))
print("DOF : ", len(body_def['modes']))

farm = project_data['farm_definition']
if farm["Ne_modes"]>0 :
    print("-- with evanescent waves (Ne modes =", farm["Ne_modes"], ")")
if farm['RAO'] :
    print("-- RAO calculation activated")

## generate wec (diffraction and force transfer matrices and radiation coefficents)
WEC = Body(freqs, directions*np.pi/180., depth, farm["Ne_modes"], convention)

if farm["Ne_modes"]>0 : 
    Evanescent=True
    WEC.Transfers(dire, results, results, Evanescent, BEM='N', Tol= 1e-6)
    WEC.PickTransfers(picklefile_E, save=True)
    loaded_WEC = WEC.PickTransfers(picklefile_E, save=False)
else : 
    # if not UT.file_exist(picklefile):
    Evanescent=False
    WEC.Transfers(dire, results, results, Evanescent, BEM='N', Tol= 1e-9)
    WEC.PickTransfers(picklefile, save=True)
    
    loaded_WEC = WEC.PickTransfers(picklefile, save=False)
WEC.Write(results, Evanescent)


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
    WECArr.Scattering(coord, Evanescent) #  WECArr.Fex.shape (Num freq, Num dir, Num dof)
    WECArr.Radiation(coord, Evanescent) # Available: WECArr.Madd and WECArr.Crad, which stand for added mass and radiation damping, respectively: WECArr.Madd.shape (Num freq, Num dof, Num dof)
    
    if farm['RAO'] :
        if not os.path.exists(Motion):
            os.makedirs(Motion)
        Mechanics = os.path.join(dire, "Mechanics")
        directoryH=Nemoh.InputHydro(dire, coord)
        statusH = Nemoh.RunHydro(directoryH, nemohdynamics)
        WECArr.RAO(Mechanics) 
    if farm['Kochin']['number']>0 :
        if not os.path.exists(Motion):
            os.makedirs(Motion)
        thetaK = np.linspace(farm['Kochin']['min'],
                        farm['Kochin']['max'],
                        farm['Kochin']['number'])
        if farm['Kochin']['format']=='DEG' :
            thetaK=thetaK*np.pi/180
        # print(thetaK)
        WECArr.Kochin(thetaK) 

    if farm['Free_surface']['Nx'] > 0 :
        if not os.path.exists(Motion):
            os.makedirs(Motion)
        WECArr.FreeSurface(coord, farm['Free_surface']['Nx'], farm['Free_surface']['Ny'], farm['Free_surface']['Lx'], farm['Free_surface']['Ly']) 


#################################################################
#################################################################
######################### OUTPUT ################################
#################################################################
#################################################################
    if not os.path.exists(resultsIT):
        os.makedirs(resultsIT)
    out=output.WriteData(WECArr, directionMB, farm, distance, resultsIT, results_h5, Motion)
    if len(farm['Plot_DOF'])>0 :
        if not os.path.exists(PlotFolder):
            os.makedirs(PlotFolder)
        out=plot.PlotData(WECArr, farm, distance, directionMB, farm['Plot_DOF'], PlotFolder)
    param_distance=param_distance*2


assert len(WECArr.Fex[:, 0, 0]) == len(WECArr.Madd[:, 0, 0]) == len(WECArr.period)
assert len(WECArr.Fex[0, 0, :]) == len(WECArr.Madd[0, :, 0]) == len(WECArr.Madd[0, 0, :]) == N_bodies*len(body_def['modes'])
assert len(WECArr.Fex[0, :, 0]) == len(directionMB)

## NOTE: Notice you only need the instance WEC to run the DIRECT MATRIX METHOD.
## Once you've got WEC you can skip the part SINGLE BODY HYDRODYNAMICS and go straight to DIRECT MATRIX METHOD for different array configurations (coord) or different wave headings (directionMB).
# plt.show()
print('end of script')

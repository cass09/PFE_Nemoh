"""
"""

# Imports
from math import pi
import numpy as np
from glob import glob
import os
import shutil
from toolbox.CalaixSastre import g, ro, execute
from toolbox.BodyMesh import get_panel_centers_cylindrical
from toolbox import mesh as _msh

# Shortcuts
_md = os.mkdir
_j = os.path.join

# Settings
IRF = [0,0.1,10] # IRF calculation (1 for yes), time step and duration
P = 0 # Save pressure on body surface (1 for yes)
Kochin = [0,0.,180.] # Kochin function calculation: Number of directions of calculation (0 for no calculations), Min and Max (degrees)
FS = [0,50,200.,200.] # Free surface visualization: Number of points in x direction (0 for no calculations) and y direction and dimensions of domain in x and y direction
RAO=0

# Script
def InputHydro(directory, coord) :
    Nb = len(coord)
    with open(_j(directory,'Hydro.txt'), 'w') as f :
        f.write('{:}				! Number of bodies\n'.format(Nb))   
        for nb in range(Nb):  
            f.write('{:} {:} {:}		! CoG\n'.format(coord[nb, 0], coord[nb, 1], 0.0))
    return directory

def InputDynamics(directory, meshes, modes, freq, directions, depth, FieldPoints, solver, sources) :
    """
    Generate input files for further running of Nemoh preProcessor.exe,
    Solver.exe and postProcessor.exe.

    directory (str): directory where to generate input files.
    meshes (list): each element is a list of: (str) mesh path for each body,
                   (int) number of nodes and (int) number of panels.
                   If given as a local path, it should be provided with respect to
                   directory !!!
    modes (list): each element is a list of: (list) each element is a (list) of
                  (int)+(6 floats) describing the generalized degrees of
                  freedom in Nemoh.cal and (list) each element is a
                  (list) of (int)+(6 floats) describing the generalized forces
                  in Nemoh.cal. modes[body][DoFs/Forces] = [(int),6x(float)].
    freq (list): [(int) Number of freqs, (float) min freq, (float) max freq].
    directions (list): [(int) Number of directions, (float) min direction,
                       (float) max direction].
    depth (float): water depth.
    FieldPoints (list): [(float) radius of the cylinder, (int) number of azimuths,
                        (int) number of axial z-coordinates]
    solver (integer) : type solver for Nemoh resolution
    """
    #
    Nb = len(modes)

    # if directory does not exist, make it
#    cannot = True
#    num = 1
#    while cannot :
#        try :
#            _md(directory + str(num))
#            cannot = False
#        except :
#            num += 1
#    directory += str(num)

    # make a mesh directory and a results directory
    _md(_j(directory, 'mesh'))
    _md(_j(directory, 'results'))

    # Generate ID.dat
    # with open(_j(directory,'ID.dat'), 'w') as f :
    #     f.write('1\n')
    #     f.write('.')

    # Generate input.txt
    with open(_j(directory,'input_solver.txt'), 'w') as f :
        f.write('2				! Gauss quadrature (GQ) surface integration, N^2 GQ Nodes, specify N(1,4)\n')
        f.write('0.001			! eps_zmin for determine minimum z of flow and source points of panel, zmin=eps_zmin*body_diameter\n')
        f.write('{:} 				! 0 GAUSS ELIM.; 1 LU DECOMP.: 2 GMRES	!Linear system solver\n'.format(solver))
        f.write('10 1e-5 1000  	! Restart parameter, Relative Tolerance, max iter -> additional input for GMRES')
    
    with open(_j(directory,'Hydro.txt'), 'w') as f :
        f.write('{:}				! Number of bodies\n'.format(Nb))   
        for nb in range(Nb):  
            f.write('0.0    0.0     0.0		! CoG\n')

    # Generate Nemoh.cal
    with open(_j(directory,'Nemoh.cal'), 'w') as f :
        f.write('--- Environment ------------------------------------------------------------------------------------------------------------------\n')
        f.write('{:.2f}				! RHO 			KG/M**3\n'.format(ro))
        f.write('{:.2f}				! G			M/S**2\n'.format(g))
        f.write('{:.2f}				! WATER DEPTH			M (0 for infinite water depth)\n'.format(abs(depth)))
        f.write('0.	0.			    ! COORDINATES OF WAVE MEASUREMENT POINT         M	\n')
        f.write('--- Description of floating bodies -----------------------------------------------------------------------------------------------\n')
        f.write('{:}				! Number of bodies\n'.format(Nb))
        for nb in range(Nb):
            f.write('--- Body {:} -------------------------------------------------------------------------------------------------------------------\n'.format(nb))
            f.write('"{:}"                ! Name of mesh file\n'.format(meshes[nb][0]))
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
        f.write('1 {:} {:} {:}		! Number of wave frequencies, min, max (rad/s)\n'.format(*freq))
        f.write('{:} {:} {:}		! Number of wave direction, min, max (degrees)\n'.format(*directions))
        f.write('--- Post processing ---------------------------------------------------------------------------------------------------------------\n')
        f.write('{:} {:} {:}		! IRF calculation (1 for yes), time step and duration\n'.format(*IRF))
        f.write('{:}				! Save pressure on body surface (1 for yes)\n'.format(P))
        f.write('{:} {:} {:}		! Kochin function calculation: nb of direction (0 for no), min, max (degrees)\n'.format(*Kochin))
        f.write('{:} {:} {:} {:}	! Free surface vizualization: nb of points in x (0 for no), y and dimensions of domain in x and y directions\n'.format(*FS))
        f.write('{:}	            ! Response Amplitude Operator (RAO), 0 no calculation, 1 calculated -> Inertia.cal\n'.format(RAO))
        f.write('1	            ! output freq type, 1,2,3=[rad/s,Hz,s]\n')
        f.write('--- Interaction Theory ---------------------------------------------------------------------------------------------------------------\n')
        f.write('{:} {:} {:}		! Field points calculation: radius of the cylinder (0 for no calculations), number of points in theta and number of points in z\n'.format(*FieldPoints))
        f.write('{:} 0		! Field points calculation: radius of the cylinder (0 for no calculations), number of points in theta and number of points in z\n'.format(sources))
        f.write('0	                ! run Interaction Theory \n')
        f.write('--- QTF ---------------------------------------------------------------------------------------------------------------\n')
        f.write('0  				! QTF flag, 1 is calculated \n')
        f.write('\n')

    return directory

def InputStatics(directory, meshes, cg, targNp, sym = 0, disp = [0.,0.]) :
    """
    Generate input files for further running of Nemoh Mesh.exe.

    directory (str): directory where to generate input files.
    meshes (list): each element is a (str) mesh path for each body.
                   If given as a local path, it should be provided with respect to
                   that used for providing the directory !!!
                   Possible formats are .gdf, .dat prior Mesh.exe, .dat and
                   SurF objects.
    cg (list): each element is a (list) with [(float) x-coordinat, (float)y-coord]
               of gravity centre.
    TargNp (list): each element is an (int) with the target for the number of
                  panels in refined mesh.
    sym (int): 1 if a symmetry about (xOz) is used. 0 otherwise.
    disp (list): Possible translation about x axis (first element (float)) and
                 y axis (second element (float)).
    """
    #
    sep = _j(' ',' ').split()[0]

    # if directory does not exist, make it
#    cannot = True
#    num = 1
#    while cannot :
#        try :
#            _md(directory + str(num))
#            cannot = False
#        except :
#            num += 1
#    directory += str(num)

    directories = list()
    for nb, mesh in enumerate(meshes) :
        # Check on mesh file
        if type(mesh) == str :
            if mesh[-3:] == 'GDF' or mesh[-3:] =='gdf' :
                Geom = _msh.readGDF(mesh)
                Geom.dat(mesh.split('.GDF')[0].split('.gdf')[0], refined = False)
                mesh = Geom.dirDAT
                status = True
            elif mesh[-3:] == 'DAT' or mesh[-3:] =='dat' :
                with open(mesh, 'r') as F :
                    line = int(F.readline().split('!')[0].split()[0])
                if line == 2 :
                    Geom = _msh.readDAT(mesh)
                    Geom.dat(mesh.split('.DAT')[0].split('.dat')[0]+'2', refined = False)
                    mesh = Geom.dirDAT
                status = True
            else :
                status = False
        elif type(mesh) == _msh.SurF :
            mesh.dat(_j(directory,'SurF{:}'.format(nb)), refined = False)
            mesh = mesh.dirDAT
            status = True
        else :
            status = False

        # Inputs
        if status == True :
            #
            name = mesh.split(sep)[-1]
            directories.append(_j(directory,'body'+str(nb+1)))

            # make a mesh directory and copy the mesh file in it
            _md(directories[nb])
            _md(_j(directories[nb],'mesh'))
            shutil.copy(mesh, _j(directories[nb],'mesh'))

            # Generate ID.dat
            with open(_j(directories[nb],'ID.dat'), 'w') as f :
                f.write('1\n')
                f.write('.')

            # Generate Mesh.cal
            with open(_j(directories[nb],'Mesh.cal'), 'w') as f :
                f.write('"{:}"		   ! name of mesh file\n'.format(name))
                f.write('{:}          ! 1 if a symmetry about x0z is used. 0 otherwise\n'.format(sym))
                f.write('{:} {:}      ! Possible translation about x axis (first number) and y axis (second number)\n'.format(*disp))
                f.write('{:} {:} {:}  ! Coordinates of gravity centre\n'.format(*cg[nb]))
                f.write('{:}          ! Target for the number of panels in refined mesh\n'.format(targNp[nb]))
                f.write('2\n')
                f.write('0.\n')
                f.write('1.\n')
                f.write('{:}          ! Water density (kg/m3)\n'.format(ro))
                f.write('{:}          ! Gravity (m/s2)\n'.format(g))

    return directories

def RunDynamics(directory, nemohdynamics = 'C:\\nemoh') :
    """
    Run Nemoh preProcessor.exe, Solver.exe and postProcessor.exe.

    directory (str): path where input files.
    nemohdynamics (str): path where .exe files.
    """
    # get current working directory
    cwd = os.getcwd()

    # change to directory where there are the inputs for dynamics
    os.chdir(directory)

    # CML : extension .exe pour Windows
    # preProcessor
    # execute(_j(nemohdynamics,'preProcessor.exe'))
    execute(_j(nemohdynamics,'preProc'))

    # execute(_j(nemohdynamics,'hydrosCal'))

    # Solver
    # execute(_j(nemohdynamics,'Solver.exe'))
    execute(_j(nemohdynamics,'solver'))

    # postProcessor
    # execute(_j(nemohdynamics,'postProcessor.exe'))
    execute(_j(nemohdynamics,'postProc'))

    # change back to initial directory
    os.chdir(cwd)

    return True

def RunHydro(directory, nemohdynamics = 'C:\\nemoh') :
    """
    Run Nemoh preProcessor.exe, Solver.exe and postProcessor.exe.

    directory (str): path where input files.
    nemohdynamics (str): path where .exe files.
    """
    # get current working directory
    cwd = os.getcwd()

    # change to directory where there are the inputs for dynamics
    os.chdir(directory)

    # CML : extension .exe pour Windows
    
    execute(_j(nemohdynamics,'hydrosCal'))
    # execute(_j(nemohdynamics,'postProc'))

    # change back to initial directory
    os.chdir(cwd)

    return True


def RunStatics(directories, nemohstatics = 'C:\\nemoh') :
    """
    Run Nemoh Mesh.exe.

    directories (list): each element is a (str) of: path where input files
                        for each body are located.
    nemohstatics (str): path where Mesh.exe.
    """
    # get current working directory
    cwd = os.getcwd()

    for directory in directories :

        # change to directory where there are the inputs for statics
        os.chdir(directory)
        # Mesh
        # execute(_j(nemohstatics,'Mesh.exe'))
        execute(_j(nemohstatics,'Mesh'))

        # change back to initial directory
        os.chdir(cwd)

    return True

def ReadDynamics(directory, DOF, FORCE, Nd, convention = 'N') :
    """
    Read Hydrodynamic Forces CA.dat, CM.dat and ExcitationForce.tec

    directory (str): path where CA.dat, CM.dat and ExcitationForce.tec.
    DOF (int): Total number of generalized degrees of freedom.
    FORCE (int): Total number of generalized forces.
    Nd (int): Total number of wave directions for ExcitationForce.tec.
    convention (str): 'N', default, for nemoh convention exp(-1j*omega*time) and
                      'W' for wamit convention exp(1j*omega*time).
    """
    # Read CA.dat
    if len(glob(_j(directory,'CA.dat'))) > 0 :
        with open(_j(directory,'CA.dat'), 'r') as f_CA :
            Nf = np.array(f_CA.readline().split(':')[-1], dtype = int)
            Ca = np.zeros((Nf, FORCE, DOF), dtype = float)
            Nlines = int(np.ceil(FORCE/6.))
            f_CA.readline()
            f_CA.readline() # CML 
            for f in range(Nf):
                f_CA.readline()
                for force in range(FORCE):
                    ls_tmp = list()
                    for nls in range(Nlines):
                        ls_tmp += f_CA.readline().split()
                    Ca[f,force,:] = np.array(ls_tmp, dtype = float)
    else : Ca = np.array([])
    # Read CM.dat
    if len(glob(_j(directory,'CM.dat'))) > 0 :
        with open(_j(directory,'CM.dat'), 'r') as f_CM :
            Cm = np.zeros((Ca.shape) , dtype = float)
            f_CM.readline()
            for f in range(Nf):
                f_CM.readline()
                for force in range(FORCE):
                    ls_tmp = list()
                    for nls in range(Nlines):
                        ls_tmp += f_CM.readline().split()
                    Cm[f,force,:] = np.array(ls_tmp,dtype = float)
    else : Cm = np.array([])
    # Read ExcitationForce.tec
    if len(glob(_j(directory,'ExcitationForce.tec'))) > 0 :
        with open(_j(directory,'ExcitationForce.tec'), 'r') as f_Ex :
            Nd = int(Nd)
            Fex = np.zeros((Nf, Nd, FORCE), dtype = complex)
            Nlines = int(np.ceil(((FORCE*2)+1)/80.))
            for skip_header in range(FORCE+1):
                f_Ex.readline()
            for d in range(Nd):
                f_Ex.readline()
                for f in range(Nf):
                    ls_tmp = list()
                    for nls in range(Nlines):
                        ls_tmp += f_Ex.readline().split()
                    MagPha = np.array(ls_tmp[1:],dtype = float)
                    Fex[f,d,:] = MagPha[0:-1:2]*np.exp(1j*MagPha[1::2])
    else : Fex = np.array([])
    if convention == 'W' :
        Fex = np.conj(Fex)
    return (Fex, Cm, Ca)

def ReadFieldPoints(directory, DOF, freq, Nd, FieldPoints, convention = 'N') :
    """
    Read Pressure from cylsurface. Num.dat

    directory (str): path where cylsurface. Num.dat.
    DOF (int): Total number of generalized degrees of freedom.
    freq (list): wave frequencies.
    Nd (int): Total number of wave directions for ExcitationForce.tec.
    FieldPoints (list): [(float) radius of the cylinder, (int) number of azimuths,
                        (int) number of axial z-coordinates]
    convention (str): 'N', default, for nemoh convention exp(-1j*omega*time) and
                      'W' for wamit convention exp(1j*omega*time).
    """
    Nf = len(freq)
    Np = DOF*Nf + Nd*Nf # Radiation + Scattering
    radius = FieldPoints[0]
    (Nth, Nz) = int(FieldPoints[1]), int(FieldPoints[2])
    Pressure = np.zeros((Nth, Nz, Np), dtype = complex)
    FieldPoints = np.zeros((Nth, Nz, 3))
    for p in range(Np):
        # print('Unpack problem {:} / {:}'.format(p+1, Np))
        with open(_j(directory,'cylsurface.{:05d}.dat'.format(p+1)), 'r') as f_prob :
            burn = f_prob.readline()
            burn = f_prob.readline()
            for nz in range(int(Nz)):
                for nth in range(int(Nth)):
                    temp = np.array(f_prob.readline().split(),dtype = float)
                    FieldPoints[nth,nz,:] = temp[0:3]
                    Pressure[nth,nz,p] = temp[3]*np.exp(1j*temp[4])
    # rearrange Pressure into Radiation or Scattering
    Scattering = np.zeros((Nf, Nd, Nz, Nth), dtype = complex)
    Radiation =  np.zeros((Nf, DOF, Nz, Nth), dtype = complex)
    ip = -1
    for f in range(Nf):
        for d in range(Nd):
            ip +=1
            Scattering[f,d,:,:] = Pressure[:,:,ip].T
        for dof in range(DOF):
            ip +=1
            Radiation[f,dof,:,:] = Pressure[:,:,ip].T
    # replace nan-values by zero
    Scattering[np.isnan(Scattering)] = 0.
    Radiation[np.isnan(Radiation)] = 0.
    # get FieldPoints
    z = FieldPoints[0,:,2]
    th = np.linspace(0, 2*pi, Nth, endpoint = False)
    del(burn)
    # PhiR_W = 1j*freq*conj(Phirad_N)
    for ind, fr in enumerate(freq):
        Scattering[ind] *= g/1j/fr
        Radiation[ind] *= g/1j/fr*(-1j*fr) # *(-1j*fr) Radiation will be due to unit amplitude motion
    if convention == 'W' :
        Scattering = np.conj(Scattering)
        Radiation = np.conj(Radiation)
    return (Scattering, Radiation, radius, th, z)

def ReadSources(directory, DOF, freq, Nd, meshfile, convention = 'N') :
    """
    Read sources terms from sources. Num.dat

    directory (str): path where cylsurface. Num.dat.
    DOF (int): Total number of generalized degrees of freedom.
    freq (list): wave frequencies.
    Nd (int): Total number of wave directions for ExcitationForce.tec.
    FieldPoints (list): [(float) radius of the cylinder, (int) number of azimuths,
                        (int) number of axial z-coordinates]
    convention (str): 'N', default, for nemoh convention exp(-1j*omega*time) and
                      'W' for wamit convention exp(1j*omega*time).
    """
    Nf = len(freq)
    Np = DOF*Nf + Nd*Nf # Radiation + Scattering
    # radius = FieldPoints[0]
    # (Nth, Nz) = int(FieldPoints[1]), int(FieldPoints[2])
    centers, areas=get_panel_centers_cylindrical(os.path.join(directory, '..', meshfile))
    Npanels=len(centers)
    Sources = np.zeros((Npanels, Np), dtype = complex)
    # FieldPoints = np.zeros((Nth, Nz, 3))
    for p in range(Np):
        with open(_j(directory,'sources.{:05d}.dat'.format(p+1)), 'r') as f_prob :
            for i in range(int(Npanels)):
                    temp = np.array(f_prob.readline().split(),dtype = float)
                    Sources[i,p] = temp[0]+1j*temp[1]
    # rearrange Pressure into Radiation or Scattering
    Scattering = np.zeros((Nf, Nd, Npanels), dtype = complex)
    Radiation =  np.zeros((Nf, DOF, Npanels), dtype = complex)
    ip = -1
    for f in range(Nf):
        for d in range(Nd):
            ip +=1
            Scattering[f,d,:] = Sources[:,ip].T
        for dof in range(DOF):
            ip +=1
            Radiation[f,dof,:] = Sources[:,ip].T
    # replace nan-values by zero
    Scattering[np.isnan(Scattering)] = 0.
    Radiation[np.isnan(Radiation)] = 0.
    # radius = centers[:,0]
    # theta = centers[:,1]
    # z = centers[:,2]
    return (Scattering, Radiation, centers, areas)

def ReadStatics(directories, whichDOFs, whichFORCEs) :
    """
    Read Hydrostatic forces KH.dat, Inertia_hull.dat and Hydrostatics.dat

    directories (list): each element is a (str) consisting of path where KH.dat,
                        Inertia_hull.dat and Hydrostatics.dat for each of the bodies.
    whichDOFs (list): each element is a (tuple) consisting of which of the usual 6
                     rigid body degrees of freedom are considered for each body.
                     The indexing should be:
                      0 "surge"
                      1 "sway"
                      2 "heave"
                      3 "roll"
                      4 "pitch"
                      5 "yaw"
    whichFORCEs (list): same as whichDOFs, but with forces instead of degrees of
                        freedom.
    """
    (KHs , Ms, rows, cols) = (list(), list(), list(), list())
    # read KH.dat
    for directory, whichDOF, whichFORCE in zip(directories, whichDOFs, whichFORCEs) :
        with open(_j(directory,'KH.dat'), 'r') as f_KH :
            KH = list()
            for line in map(f_KH.readlines().__getitem__ , whichFORCE) :
                KH.append(np.array(line.split() , dtype = float)[whichDOF,])
        KHs.append(np.array(KH, dtype = float))
        # read Inertia_hull.dat
        with open(_j(directory,'Inertia_hull.dat'), 'r') as f_I :
            I = list()
            for line in f_I.readlines() :
                I.append(np.array(line.split() , dtype = float))
        I = np.array(I, dtype = float)
        # read Hydrostatics.dat
        with open(_j(directory,'Hydrostatics.dat'), 'r') as f_H :
            Vol = float(f_H.readlines()[3].split('=')[-1])
        # mass matrix
        M = np.zeros((6,6), dtype = float)
        M[(0,1,2),(0,1,2)] = Vol*ro
        M[3:,3:] = I
        Ms.append(M[whichFORCE,:][:,whichDOF])
        # number of rows
        rows.append(len(whichFORCE))
        cols.append(len(whichDOF))
    KH = np.zeros((sum(rows),sum(cols)), dtype = float)
    M = np.zeros((sum(rows),sum(cols)), dtype = float)
    for nb in range(len(directories)) :
        rs = [sum(rows[:nb]) , sum(rows[:nb]) + rows[nb]]
        cs = [sum(cols[:nb]) , sum(cols[:nb]) + cols[nb]]
        KH[rs[0]:rs[1],cs[0]:cs[1]] = KHs[nb]
        M[rs[0]:rs[1],cs[0]:cs[1]] = Ms[nb]
    return (KH , M)

def ReadNemohcal(directory) :
    """
    Read setting file Nemoh.cal.

    directory (str): path where Nemoh.cal or nemoh.cal
    """
    # dictionary
    dictionary = np.array([[1,1,0,0],
                           [1,0,1,0],
                           [1,0,0,1],
                           [2,1,0,0],
                           [2,0,1,0],
                           [2,0,0,1]],dtype = float)
    # read Nemoh.cal
    whichDOF = list()
    whichFORCE = list()
    DOF = 0
    FORCE = 0
    with open(_j(directory , 'Nemoh.cal') , 'r') as f_N :
        for il , line in enumerate(f_N) :
            if il ==  3 :
                WaterDepth = float(line.split('!')[0])
            elif il == 6 :
                Nb = int(line.split('!')[0]) # Number of bodies
            elif il > 6 :
                if il <= 6 + Nb :
                    whichDOFb = list()
                    for ilbd , linebd in enumerate(f_N) :
                        if ilbd == 0 :
                            # BodyMesh = np.array(linebd.split('!')[0].split(), dtype = float)
                            BodyMesh = linebd.split('!')[0].strip().strip('"')
                        if ilbd == 2 :
                            DOFb = int(linebd.split('!')[0])
                            DOF += DOFb
                        if ilbd > 2 :
                            if ilbd <= 2 + DOFb :
                                dof = np.array(linebd.split('!')[0].split(), dtype = float)
                                for di , d in enumerate(dictionary) :
                                    if all(d == dof[:-3]) :
                                        whichDOFb.append(di)
                            if ilbd > 2 + DOFb :
                                whichDOF.append(whichDOFb)
                                FORCEb = int(linebd.split('!')[0])
                                FORCE += FORCEb
                                break
                    whichFORCEb = list()
                    for ilbf , linebf in enumerate(f_N) :
                        if ilbf < FORCEb :
                            force = np.array(linebf.split('!')[0].split(), dtype = float)
                            for fi , f in enumerate(dictionary) :
                                if all(f == force[:-3]) :
                                    whichFORCEb.append(fi)
                        if ilbf == FORCEb :
                            whichFORCE.append(whichFORCEb)
                            extralines = int(linebf.split('!')[0])
                            if extralines > 0 :
                                for iex , extra in enumerate(f_N) :
                                    if iex == extralines - 1 :
                                        break
                            break
                elif il > 6 + Nb :
                    if  il == 6 + Nb + 2 :
                        (type, Nf,fmin,fmax) = np.array(line.split('!')[0].split(), dtype = float)
                        freq = np.linspace(fmin, fmax, int(Nf))
                    elif il == 6 + Nb + 3 :
                        (Nd,dmin,dmax) = np.array(line.split('!')[0].split(), dtype = float)
                        directions = np.linspace(dmin, dmax, int(Nd))*pi/180
                    elif il == 6 + Nb + 12 :
                        FieldPoints = np.array(line.split('!')[0].split(), dtype = float)
                    elif il == 6 + Nb + 13 :
                        Method = np.array(line.split('!')[0].split(), dtype = float)
                        break
    return (DOF, FORCE, freq, directions, FieldPoints, BodyMesh, whichDOF, whichFORCE, WaterDepth, Method)

def ReadNemohcal_mshs_mds(directory) :
    """
    Read setting file Nemoh.cal and return mesh files and generalized modes.

    directory (str): path where Nemoh.cal or nemoh.cal
    """
    meshes = list()
    modes = list()
    with open(os.path.join(directory , 'Nemoh.cal') , 'r') as f_N :
        content = f_N.readlines()[6:]
        Nb = int(content[0].split('!')[0]) # Number of bodies
        strt = 1
        for body in range(Nb):
            strt += 1 #--- Body body -------------------------------------------------------------------------------------------------------------------
            mesh_file = content[strt].split('!')[0].split()[0]
            mesh_file = mesh_file.replace('"', '')
            mesh_file = mesh_file.replace("'", '')
            if mesh_file[0] == '.':
                meshes.append(directory+mesh_file[1:])
            elif meshes[0][0] == os.path.join('a','a')[1]:
                meshes.append(directory+mesh_file)
            else:
                meshes.append(os.path.join(directory, mesh_file))
            Ndof = int(content[strt+2].split('!')[0])
            DoF = [np.array(content[strt+3+ndof].split('!')[0].split(), float) for ndof in range(Ndof)]
            Nforce = int(content[strt+3+Ndof].split('!')[0])
            FoRcE = [np.array(content[strt+4+Ndof+nforce].split('!')[0].split(), float) for nforce in range(Nforce)]
            modes.append([DoF, FoRcE])
            Nadd = int(content[strt+4+Ndof+Nforce].split('!')[0])
            strt += 5+Ndof+Nforce+Nadd
    return (meshes, modes)
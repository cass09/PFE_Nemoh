"""
Read WAMIT output files
"""

from copy import copy
import numpy as np
from math import pi
from glob import glob
from toolbox.CalaixSastre import len2, WNumber, ro, g

def ReadCFGandPOT(directory) :
    """
    """
    ## Read .cfg and check IPERIN
    perORfreq = 1
    dof_gener = 0
    with open(glob(''.join((directory,'\*.cfg')))[0],'r') as fcfg:
        for line in fcfg:
            assignment = line.partition('=')
            if assignment[0].split()[0] == 'IPERIN':
                perORfreq = int(assignment[2])
            elif assignment[0].split()[0] == 'NEWMDS':
                dof_gener = int(assignment[2])
           
    ## Read .pot and get number of wave directions (Ndir) and wave periods (Nper)
    fpot = open(glob(''.join((directory,'\*.pot')))[0],'r').readlines() # pot file from WAMIT
    wdepth = float(fpot[1].split()[0])
    
    # Get periods
    Nper= int(fpot[3].split()[0]) # Nper is found in the 4th line (3 for python) and 1st column (0 for python). split() is used for breaking characters into a list with sep=whitespaces
    rowORcolarray = np.array(fpot[3+1].split(),dtype=float) # used to know whether they are arranged through a row (one line) or through column
    adjust = 1 # used to know where Ndir starts from
    if perORfreq == 2:
        if len(rowORcolarray) < Nper:
            period = 2*pi/np.array(fpot[3+1:3+1+Nper],dtype=float)
        elif len(rowORcolarray) == Nper:
            period = 2*pi/rowORcolarray
            adjust = 1./Nper
        else :
            raise IOError('Check on .POT file wave periods')
    else: 
        if len(rowORcolarray) < Nper:
            period = np.array(fpot[3+1:3+1+Nper],dtype=float)
        elif len(rowORcolarray) == Nper:
            period = rowORcolarray
            adjust = 1./Nper
        else :
            raise IOError('Check on .POT file wave periods')
            
    # Get directions
    start = 3+1+int(Nper*adjust)
    Ndir= int(fpot[start].split()[0])
    rowORcolarray = np.array(fpot[start+1].split(),dtype=float)
    adjust = 1
    if len(rowORcolarray) < Ndir:
            direction = np.array(fpot[start+1:start+1+Ndir],dtype=float)
    elif len(rowORcolarray) == Ndir:
        direction = rowORcolarray
        adjust = 1./Ndir
    else :
        raise IOError('Check on .POT file wave directions')
        
    # Get degrees of freedom info
    start += 1+int(Ndir*adjust)
    NBODY = int(fpot[start].split()[0])
    dofi_solid = np.zeros((NBODY,6), dtype=int)
    dof_solid = 0
    for body in range(NBODY):
        dofi_solid[body]= np.array(fpot[start+3+3*body].split(), dtype=int)
        dof_solid += sum(dofi_solid[body]) # Total number of degrees of freedom corresponding to regular rigid body motion (max 6)
    # dof_glob contains for each row (body) the global indexing 6 solid rigid dof + generalized dof
    dof_glob = np.array(range(1,6*NBODY+1),dtype=int).reshape((NBODY,6))
    dof_glob += np.repeat(np.array(range(NBODY))*dof_gener,6).reshape(dof_glob.shape)
    # from dof_glob I just pick the solid rigid dof accounted for each body
    dofi_solid= dof_glob[np.array(dofi_solid,dtype=bool)] 
    dof = dof_solid+dof_gener*NBODY
    del(fpot)
    return (period, direction, wdepth, NBODY, dof, dof_gener, dof_solid, dofi_solid)
            
def ReadDynamics(directory, period, direction, dof, convention = 'W'): 
    """
    """
    Nper = len(period)
    Ndir = len(direction)
    ## Read .1 (Madd and Crad)
    if len(glob(''.join((directory,'\*.1')))) > 0 :
        with open(glob(''.join((directory,'\*.1')))[0],'r') as f1: # 1 file from WAMIT
            Madd= list()
            Crad= list()
            for i,line in enumerate(f1.readlines()[1:]):# skip first line, i.e. header
                linep= line.split() 
                Madd.append(float(linep[-2]))
                Crad.append(float(linep[-1])) # Cradij component is given in the last column
        dofcheck = int(np.sqrt((i+1)/Nper))
        if dofcheck != dof:
            raise IOError('Total number of dof from Madd does not match total number of dof from dof_solid + dof_gener*NBODY')
        Madd= ro*np.array(Madd, dtype=float).reshape((Nper,dof,dof))
        Crad= np.array(Crad, dtype=float).reshape((Nper,dof,dof)) # The reshape is easy the way Crad is given by WAMIT (as it was seen for the mass matrix)
        periodp = period.repeat(dof**2).reshape(Crad.shape)
        Crad= ro*2*pi/periodp*Crad
    else : (Madd, Crad) = (np.array([]), np.array([]))

    # Read .2 (Fex)
    if len(glob(''.join((directory,'\*.2')))) > 0 :
        with open(glob(''.join((directory,'\*.2')))[0],'r') as f2: # 2 file from WAMIT
            Fex= np.zeros(Nper*Ndir*dof,dtype=complex)
            for i,line in enumerate(f2.readlines()[1:]):# skip first line, i.e. header                
                linep= line.split()
                Fex[i]= complex(float(linep[-2]),float(linep[-1])) # Fex is given as real(fex) and im(fex) in the last two columns
        Fex= ro*g*Fex.reshape((Nper,Ndir,dof)) # The reshape is easy the way Fex is given by WAMIT (as it was seen for the mass matrix)
    else : Fex = np.array([])
    if convention == 'N' :
        Fex = np.conj(Fex)
    return (Fex, Madd, Crad)
    
def ReadFieldPoints(directory, period, direction, wdepth, dof, dof_solid, convention = 'W') :
    """
    """
    Nper = len(period)
    Ndir = len(direction)
    ## Read .fpt and create fpt array and calculate r,t,z from the regarded cylinder
    with open(glob(''.join((directory,'\*.fpt')))[0],'r') as ffpt: # fpt file from WAMIT
        fpt= [] 
        for line in ffpt.readlines()[1:]:# skip first line, i.e. header 
            linep= line.split() # split characters into a list where sep is whitespace
            #linep= array(linep,dtype=float)
            fpt.append(np.array(linep,dtype=float))#[linep[0],linep[1],linep[2]),linep[3]])
    fpt= np.array(fpt) # fpt= [ID, x(ID), y(ID), z(ID)]
    r=np.sqrt(fpt[0,1]**2+fpt[0,2]**2) # radius of the cylinder
    x= fpt[fpt[:,3]==fpt[0,3],1]
    y= fpt[fpt[:,3]==fpt[0,3],2]
    t=np.arctan2(y,x) # theta of the cylinder
    z=fpt[fpt[:,1]==fpt[0,1],3] # z of the cylinder

    # Total number of field points
    Npoints= len2(z)*len2(t)
    
    ## Read directory.6p and get the number of characters within radiation and diffraction lines
    with open(glob(''.join((directory,'\*.6p')))[0],'r') as f6p: # 6p file from WAMIT
        burnheader= f6p.readline(); del(burnheader)
        lineR= len2(f6p.readline()) # get the total number of characters within the 1st line
        if dof > dof_solid:
            lineR+= len2(f6p.readline()) # if generalized modes are regarded a 2nd line should be read
        inbetween=f6p.read(lineR*(Npoints-1)); del(inbetween)
        lineD= len2(f6p.readline())
        
    ## Read directory.6p and create PhiD and PhiR
    with open(glob(''.join((directory,'\*.6p')))[0],'r') as f6p: # 6p file from WAMIT
        burnheader= f6p.readline(); del(burnheader)
        PhiD= np.zeros(Nper*Ndir*len2(z)*len2(t),dtype=complex)
        PhiR= np.zeros((Nper*dof,len2(z)*len2(t)),dtype=complex)
        for i in range(Nper):
            
            print('{:.2f} %'.format(float(i)/(Nper-1)*100))
            # Radiation
            
            linep= np.array(f6p.read(lineR*Npoints).split(),dtype=float)\
            .reshape((Npoints,2+dof*2))
            
            wfreq= 2*pi/period[i]

            dimR= g/wfreq**2
            
            PhiR[dof*i:dof*(i+1),:]= dimR*np.transpose(linep[:,2:][:,range(0,dof*2,2)]+\
            1j*linep[:,3:][:,range(0,dof*2,2)])
            
            dimD= 1j*g/wfreq
            
            # Diffraction
            for i2 in range(Ndir):
                
                linep= np.array(f6p.read(lineD*Npoints).split(),dtype=float)\
                .reshape((Npoints,5))
                
                PhiD[Npoints*Ndir*i:Npoints*Ndir*(i+1)][Npoints*i2:Npoints*(i2+1)]= dimD*(\
                linep[:,-2]+1j*linep[:,-1])
                
    PhiD= np.reshape(PhiD,(Nper,Ndir,len2(z),len2(t))) # the way the nodes of the cylinder are given and the way WAMIT provides PhiD facilitates the reshape
    PhiR= np.reshape(PhiR,(Nper,dof,len2(z),len2(t)))
    
    # Radiation due to unit amplitude motion
    for f,fr in enumerate(2*pi/period) :
        PhiR[f] *= 1j*fr
    
    direction *= pi/180

    wfreq,wdir,fz,fx= np.meshgrid(2*pi/period,direction,z,x,indexing='ij',
                               sparse=True)
    
    # check
    if abs(wdepth - max(abs(z))) > 1e-2:
        raise IOError('abs(wdepth - max(abs(z))) > 1e-2')
    wnumber= WNumber(period,wdepth)
    
    wnumber,wdir,fz,fy= np.meshgrid(wnumber,direction,z,y,indexing='ij',
                                 sparse=True)
    
    PhiP= 1j*g/wfreq*np.cosh(wnumber*(fz+wdepth))/np.cosh(wnumber*wdepth)*\
    np.exp(-1j*wnumber*(fx*np.cos(wdir)+fy*np.sin(wdir)))
    
    PhiS = PhiD-PhiP
    
    if convention == 'N' :
        PhiS = np.conj(PhiS)
        PhiR = np.conj(PhiR)    
    return (PhiS, PhiR, r, t, z)
    
    
def ReadStatics(directory, NBODY, dof, dof_gener, dofi_solid) :    
    ## Read .mmx and create M (Mass matrix)
    fmmx= open(glob(''.join((directory,'\*.mmx')))[0],'r').readlines()[6:] # mmx file from WAMIT
    Mbody= []
    M = np.zeros((dof,dof), dtype = float)
    fini = 0
    for body in range(NBODY):
        start = 44*body+3+4
        for line in fmmx[start:start+(6+dof_gener)**2]: # Mass matrix is started given in line 9
            linep= line.split()
            I = int(linep[0])
            J = int(linep[1])
            cond1 = ((6+dof_gener)*body<dofi_solid)
            cond2 = (dofi_solid<(6+dof_gener)*(body+1))
            dofibody = dofi_solid[cond1*cond2]-(6+dof_gener)*body
            if (any(I == dofibody) or I>6) and \
            (any(J == dofibody) or J>6):
                Mbody.append(linep[2]) # M[row_i,column_j]. M is given [row1,...,rowdof] so then afterwards a reshape will work super!
        ini = copy(fini)
        fini += len(dofibody)+dof_gener
        M[ini:fini,ini:fini] = np.array(Mbody,dtype=float)\
        .reshape((len(dofibody)+dof_gener,len(dofibody)+dof_gener))
        Mbody = []
    del(fmmx)


    ## Read .hst and create Khyd (Hydrostatic stiffness matrix)
    with open(glob(''.join((directory,'\*.hst')))[0],'r') as fhst: # hst file from WAMIT
        Khyd= list()
        aux = [[i1+i2 for i2 in range(dof_gener)] for i1 in range(7,6*(NBODY+1),6)]
        dofi_gener = np.array(aux, dtype=int).reshape(-1, order = 'F')            
        for line in fhst.readlines()[1:]:
            linep= line.split()
            I = int(linep[0])
            J = int(linep[1])
            if (any(I == dofi_solid) or any(I == dofi_gener)) and \
            (any(J == dofi_solid) or any(J == dofi_gener)):
                Khyd.append(linep[2])
        Khyd= ro*g*np.array(Khyd,dtype=float).reshape((dof,dof))
        
    return (Khyd, M)
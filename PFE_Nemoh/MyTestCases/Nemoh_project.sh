
distance=16
configuration=False
config_type="X"
BEM_Nb=1						# Number of bodies
Nom_projet="BEM_Hydro_16"
FS=0

# --- Description of floating bodies -  ----------------------------------------------------------------------------------------------
CdG1=(28.72 0.0 0.0)
mesh_file1="barge.dat"			# Name of mesh file
Rcyl=6.36

if [[ "$configuration" == "True" ]]; then
    source tools_bash/config.sh
else
    translate1=(0.0 0.0 0.0)
    translate2=(10 0.0 0.0)
    Dossier_Project=$Nom_projet
    lenghtX=0
    lenghtY=30
    CdG2=(0.0 0.0 0.0)
    CdG3=(0.0 0.0 0.0)
    mesh_file2="barge.dat"			# Name of mesh file
    mesh_file3="barge.dat"	
fi

N_DOF=6			# Number of degrees of freedom
N_Forces=6						# Number of resulting generalised forces
# 1 activated
DOF=(1 1 1 1 1 1)   # Surge, Sway, Haeve, Roll, Pitch, Yaw
FORCES=(1 1 1 1 1 1) # Fx, Fy, Fz, Mx, My, Mz

#  --- Load cases to be solved -------------------------------------------------------------------------------------------------------
w_type=1   
Nw=20
w_min=0.3	
w_max=2.0			# Freq type 1,2,3=[rad/s,Hz,s], Number of wave frequencies/periods, Min, and Max
BEM_Nbeta=1	
BEM_BetaMin=0.0
BEM_BetaMax=0.0				# Number of wave directions, Min and Max (degrees)

# --- Post processing ---------------------------------------------------------------------------------------------------------------
IRF=0
dt=0.1
tf=10.0				# IRF calculation (0 for no calculation), time step and duration
show_pressure=0						# Show pressure
KochinN=0
KochinMin=0.0
KochinMax=180.0			# Kochin function 		# Number of directions of calculation (0 for no calculations), Min and Max (degrees)

# Free surface elevation 	# Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction
FS_Nx=$lenghtX
FS_Ny=$lenghtY
FS_Lx=$lenghtX
FS_Ly=$lenghtY
	
RAO=1		# Response Amplitude Operator (RAO), 0 no calculation, 1 calculated -> Inertia.cal
w_type_output=1						# output freq type, 1,2,3=[rad/s,Hz,s]


# --- Environment ------------------------------------------------------------------------------------------------------------------
RHO=1000.0					# RHO 		# KG/M**3 	# Fluid specific volume 
G=9.81					# G			# M/S**2	# Gravity
DEPTH=20.0						# DEPTH			# M		# Water depth
XEFF=0.0	
YEFF=0.0					# XEFF YEFF		# M		# Wave measurement point

# ---QTF---
QTF=0         				# QTF flag, 1 is calculated 

# ---Solver----
Gauss_N=2				# Gauss quadrature (GQ) surface integration, N^2 GQ Nodes, specify N(1,4)
eps_zmin=0.001			# eps_zmin for determine minimum z of flow and source points of panel, zmin=eps_zmin*body_diameter
Solver_type=1 				# 0 GAUSS ELIM.; 1 LU DECOMP.: 2 GMRES	#Linear system solver
Restart=10 
Tol=1e-5 
MaxIter=1000  	# Restart parameter, Relative Tolerance, max iter -> additional input for GMRES

# ---Interaction Theory---
cylR=0
cylNtheta=31
cylNz=29		# Interaction Theory Cylindrical Envelop
run_IT=0			# run IT
run_BEM=1			# run BEM
IT_Nb=2			# Nb bodies
Bcoord1=(0.0   0.0)			# Coord body 1
Bcoord2=(30.0  0.0)			# Coord body 2
IT_Nbeta=1
IT_BetaMin=0.0
IT_BetaMax=0.0		# wave directions Nb, min, max

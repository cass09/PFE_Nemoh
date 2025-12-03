% Written by Sergej Antonello Sirigu 27/12/2022

% Nemoh version 3.0 ( https://gitlab.com/lheea/Nemoh )

% Inputs Nemoh:
% Fluid specific volume (RHO)
% Gravity (G)
% Water depth (DEPTH)
% Wave measurement point
% Number of bodies
% Number of points and number of panels for each body (nx(c) nf(c))
% Number of degrees of freedom
% Surge,Sway,Heave,Roll,Pitch and Yaw
% Force in x,y,z direction and Moment force in x,y,z direction
% Number of lines of additional information
% Number of wave frequencies, Min, and Max (rad/s)
% Number of wave directions, Min and Max (degrees)
% IRF calculation (0 for no calculation), time step and duration
% Show pressure
% Kochin function 		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)
% Free surface elevation 	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction
% CG(nBodies,3)     : position of gravity centre

% Outputs :
% - A  : Matrix (6xnBodies)x(6xnBodies)xlength(w) of added mass coefficients
% - B  : Matrix (6xnBodies)x(6xnBodies)xlength(w) of radiation damping coefficients
% - Fe : Matrix (6xnBodies)xlength(w) of excitation forces (complex
% values)

clc
clear all
close all
clear py;                        % Nettoyage général
py.importlib.invalidate_caches(); % Invalide les caches Python

pyenv('ExecutionMode','OutOfProcess');

%% CHANGE TO SCRIPT DIR
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

%% FLAGS
nemohInput.ID_PLOT_RESULTS = 1;
ID_QTF = 0;                        % Flag to activate QTF computation (0 or 1)
nemohInput.ID_HydrosCal=0;         % A switch,1 computes hydrostatics, inertia, kH.  % It has been modified. Another routine is used to calculate hydrostatics
nemohInput.QTFInput=[ID_QTF,2];    %[Switch,Contrib]

%% PATH DEFINITION 
% Determine where your m-file's folder is.
[pathstr,~,~] = fileparts(mfilename('fullpath'));
% Add that folder plus all subfolders to the path.
addpath(genpath(pathstr)); % Include the subfolders in the Matlab PATH

% Vérifie que Nemoh est disponible et ajoute-le au PATH pour cette session
assert(FindingNemoh(ID_QTF, true))

# TO ADAPT !!!!!!!!!!!!
insert(py.sys.path, int32(0), '/home/cassandra/.local/lib/python3.10/site-packages'); # gmsh
insert(py.sys.path, int32(0), '/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/Matlab/BEM_matlab/python_mesh');
modulesPY = {'pontoon', 'gmsh_to_nemoh_N'};
for k = 1:length(modulesPY)
    mod = py.importlib.import_module(modulesPY{k});
    py.importlib.reload(mod);
end
%% PROJECT DIRECTORY AND FILE NAME
jsonText = fileread('Nemoh_BEM_matlab.json');
Jf = jsondecode(jsonText);
projdir = sprintf('%s%s', Jf.project.folder, Jf.project.name);   % working project folder
savefilename = 'hydrodynamicResults'; 

% Folder and file name mesh input
folderNemoh = 'nemohMeshFolder';             
nameNemohMesh=Jf.mesh.filename;  

launch_Nemoh=Jf.project.run_BEM;

if Jf.mesh.create 
    if Jf.mesh.Wings(1)==0 
        py.pontoon.create_Pontoon(Jf.mesh.Pontoon,Jf.mesh.deltaX, fullfile(folderNemoh, nameNemohMesh));
    else 
        py.PontoonWings.create_PontoonWings(Jf.mesh.Pontoon,Jf.mesh.Wings,Jf.mesh.deltaX, fullfile(folderNemoh, nameNemohMesh));
    end
    py.gmsh_to_nemoh_N.convert(fullfile(folderNemoh, nameNemohMesh));
end 
%% NEMOH INPUTS

% Bodies
nBodies = Jf.farm.N_bodies;     %number of bodies
CG_orig = Jf.params.rotation_center;  % Center of gravity of bodies in Local Reference Frame (z axis water level)
nemohInput.transMatrix = zeros(nBodies, 3);
for i = 1 : nBodies
    pos = Jf.farm.layout(i).position;  % [x y]
    nemohInput.transMatrix(i,:) = [pos(1), pos(2), 0];  % z = 0
end                                  
% nemohInput.transMatrix = [0 0 0];  % Translation matrix
nemohInput.CG = CG_orig.' + nemohInput.transMatrix;  
% Environment
nemohInput.rho = 1000;      % water density
nemohInput.g = 9.81;     %gravity
nemohInput.depth = Jf.params.water_depth;      %sea depth% Center of gravity of bodies il Global Reference Frame

% Nemoh setup 
nemohInput.wavetype=Jf.params.frequency.format;   % Freq type 1,2,3=[rad/s,Hz,s]
nemohInput.wmax=Jf.params.frequency.max;     % max frequency (rad/s)
nemohInput.wmin = Jf.params.frequency.min;  % min frequency (rad/s)
nemohInput.nbfreq=Jf.params.frequency.number;       % number frequencies
nemohInput.wavefreq= linspace(nemohInput.wmin,nemohInput.wmax,nemohInput.nbfreq)';  % min w [rad/s], max w [rad/s], Nw
nemohInput.wavedir = [Jf.params.direction.number	Jf.params.direction.min	Jf.params.direction.max];                    %Number of wave directions, Min and Max (degrees)
nemohInput.IRFcalculation = [0	0.05	60.];              %IRF calculation (0 for no calculation), time step and duration
nemohInput.showpressure = 0;                              %Show pressure
nemohInput.kochinfunction = [0	0.	180.];                 %Kochin function 		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)
nemohInput.freesurfaceelevation = [Jf.params.Free_surface.Nx	Jf.params.Free_surface.Ny	Jf.params.Free_surface.Lx	Jf.params.Free_surface.Ly]; %Free surface elevation 	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction
nemohInput.DOF = Jf.params.DOF;
nemohInput.nDegFreedom = sum(nemohInput.DOF);
nemohInput.nBodies = nBodies;
nemohInput.RAO = Jf.params.RAO;
nemohInput.OuterCylinder = [Jf.params.cylinder.radius  10  10];
nemohInput.Source = Jf.params.sources;

%% Replicate and translate body mesh (can be modified for different bodies)
for ii = 1:nBodies
nameNemohMesh_trasl{ii,1} =  [nameNemohMesh,'_',num2str(ii)];  % name of current .dat nemoh mesh 
[nodes(ii,1), panels(ii,1)] = translateMesh(nBodies,folderNemoh,nameNemohMesh,nameNemohMesh_trasl{ii,1},nemohInput.transMatrix(ii,:)); 
end

nemohInput.nodes = nodes;
nemohInput.panels = panels;
%% Create Nemoh project directory

if exist(projdir, 'dir')
      rmdir(projdir,'s') 
end
mkdir(projdir);
%system(['mkdir ',projdir]);
system(['chmod u+rwx ', projdir]); 
system(['mkdir ',projdir,filesep,'mesh']);
system(['mkdir ',projdir,filesep,'results']);

%% Moving the nemoh mesh.dat file to nemoh project folder
 for ii=1:length(nameNemohMesh_trasl)
     movefile([folderNemoh,filesep,nameNemohMesh_trasl{ii},'.dat'],[projdir],'f');
 end



%% NEMOH SOLVER PARAMETERS
writeInputSolverFile(projdir); % Nemoh Solver Parameters to edit
%% HYDROSTATICS 
% Compute hydrostatic parameters for each body 
writeHydroFile(projdir,nemohInput.nBodies,nemohInput.CG);
[KH,Inertia,XB,YB,ZB,WPA] = hydrostatic_cal(projdir,nameNemohMesh_trasl,nemohInput);
pause(1)
%% BEM NEMOH
% Write Nemoh files for simulation
% TO BE IMPLEMENTED: IRREGULAR FREQUENCIES REMOVAL
writeNemohCalFile(projdir,nameNemohMesh_trasl,nemohInput);
waveNumber(projdir, nemohInput);
if launch_Nemoh
[Idw,w,A,B,Fe] = BEM_Nemoh(projdir,nemohInput.ID_HydrosCal,ID_QTF,nBodies);
 
if nemohInput.nDegFreedom==6
    %% BUILD MATRICES
    M = zeros(6*nBodies,6*nBodies);
    K = zeros(6*nBodies,6*nBodies);
    
    % INERTIA MATRIX IS HERE CALCULATED BY NEMOH (INERTIA of the submerged displaced water volume)
    % If the mass matrix of the body is known replace this part of the code
    for ii = 1:nBodies
        M(1+(ii-1)*6:ii*6,1+(ii-1)*6:ii*6) = squeeze(Inertia(ii,:,:)); % Inertia multi-body matrix
        K(1+(ii-1)*6:ii*6,1+(ii-1)*6:ii*6) = squeeze(KH(ii,:,:));      % Stiffness multi-body matrix
    end
    
    
    RAO = zeros(length(w),6*nBodies);
    for jj = 1:length(w)
        wval = w(jj);
        RAO(jj,:) = Fe(jj,:)/(-wval^2*(M+A(:,:,jj))+1i*wval*B(:,:,jj)+K);
    end
    % save(savefilename,'M','K','A','B','Fe','w')
end

end

%% PLOT EXAMPLE TO EDIT CASE BY CASE
% DoF = 5
% figure
% plot(w,abs(RAO(:,DoF))*180/pi,'-sk')
% hold on
% DoF = DoF+6
% plot(w,abs(RAO(:,DoF))*180/pi,'-dr')
% DoF = DoF+6
% plot(w,abs(RAO(:,DoF))*180/pi,'-xb')

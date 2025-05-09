%--------------------------------------------------------------------------------------
%
%    Copyright (C) 2022 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%    Contributors: A. Babarit, R. Kurnia
%-------------------------------------------------------------------------
%
% --> function axiMesh(r,z,n,nang,npanelt,CoG,depth,wavefreq,wavedir,
%                       QTFInput,projdir)
%
% Purpose : Mesh generation of an axisymmetric body for use with Nemoh
%
% Inputs : description of radial profile of the body
%   - r         : array of radial coordinates
%   - z         : array of vertical coordinates
%   - n         : number of radial points for discretisation
%   - nang      : number of radial points for discretisation
%   - npanelt   : number of panels target
%   - CoG       : position of centre of gravity, [x,y,z]
%   - depth     : mean water depth
%   - wavefreq  : wave frequency [rad/s] input, [w1,w2,w3,...]
%   - wavedir   : wave direction [deg] input, [beta1,beta2,beta3,...]
%   - QTFInput  : QTF parameters input; [Switch 0 or 1, Contrib]
%   - projdir   : path of the project directory
%
% Outputs : produces NEMOH input files in the project directory. The files
%           are Mesh.cal, Nemoh.cal,input_solver.txt,[meshfile].dat,
%           hydrostatic stiffness matrix KH.dat,Inertia.dat,Hydrostatics.dat
%
% Warning : z(i) must be greater than z(i+1)
%
%
function axiMesh(r,z,n,nang,npanelt,CoG,depth,wavefreq,wavedir,QTFInput,projdir)
status=close('all');
system(['mkdir ',projdir]);
system(['mkdir ',projdir,filesep,'mesh']);
system(['mkdir ',projdir,filesep,'results']);
nBodies=1;
ntheta=nang;%input('\n - Enter number of points for angular discretisation : ');
theta=[0.:pi/(ntheta-1):pi];
nx=0;
% Calcul des sommets du maillage
for j=1:ntheta
    for i=1:n
        nx=nx+1;
        x(nx)=r(i)*cos(theta(j));
        y(nx)=r(i)*sin(theta(j));
        z(nx)=z(i);
    end
end
% Calcul des facettes
nf=0;
for i=1:n-1
    for j=1:ntheta-1
        nf=nf+1;
        NN(1,nf)=i+n*(j-1);
        NN(2,nf)=i+1+n*(j-1);
        NN(3,nf)=i+1+n*j;
        NN(4,nf)=i+n*j;
    end
end
% Affichage de la description du maillage
nftri=0;
for i=1:nf
    nftri=nftri+1;
    tri(nftri,:)=[NN(1,i) NN(2,i) NN(3,i)];
    nftri=nftri+1;
    tri(nftri,:)=[NN(1,i) NN(3,i) NN(4,i)];
end
figure;
trimesh(tri,x,y,z,[zeros(nx,1)]);
title('Characteristics of the discretisation');
fprintf('\n --> Number of nodes             : %g',nx);
fprintf('\n --> Number of panels (max 2000) : %g \n',nf);
% Creation des fichiers de calcul du maillage
fid=fopen([projdir,filesep,'Mesh.cal'],'w');  %
fprintf(fid,'axisym \n',1);
fprintf(fid,'1 \n 0. 0. \n ');
zG=CoG(3);%input(' - Vertical position of gravity center : ');
fprintf(fid,'%f %f %f \n',[0. 0. zG]);
nfobj=npanelt;%input(' - Target for number of panels : ');
fprintf(fid,'%g \n 2 \n 0. \n 1.\n',nfobj);
fprintf(fid,'%f \n',1025.);
fprintf(fid,'%f \n',9.81);
status=fclose(fid);

fid=fopen([projdir,filesep,'axisym'],'w');
fprintf(fid,'%g \n',nx);
fprintf(fid,'%g \n',nf);
for i=1:nx
    fprintf(fid,'%E %E %E \n',[x(i) y(i) z(i)]);
end
for i=1:nf
    fprintf(fid,'%g %g %g %g \n',NN(:,i)');
end
status=fclose(fid);
% Raffinement automatique du maillage et calculs hydrostatiques
system(['mesh ',projdir]);
% Visualisation du maillage
clear x y z NN nx nf nftri tri u v w;
fid=fopen([projdir,filesep,'mesh',filesep,'axisym.tec'],'r');
ligne=fscanf(fid,'%s',2);
nx=fscanf(fid,'%g',1);
ligne=fscanf(fid,'%s',2);
nf=fscanf(fid,'%g',1);
ligne=fgetl(fid);
fprintf('\n Characteristics of the mesh for Nemoh \n');
fprintf('\n --> Number of nodes : %g',nx);
fprintf('\n --> Number of panels : %g\n \n',nf);
for i=1:nx
    ligne=fscanf(fid,'%f',6);
    x(i)=ligne(1);
    y(i)=ligne(2);
    z(i)=ligne(3);
end
for i=1:nf
    ligne=fscanf(fid,'%g',4);
    NN(1,i)=ligne(1);
    NN(2,i)=ligne(2);
    NN(3,i)=ligne(3);
    NN(4,i)=ligne(4);
end
nftri=0;
for i=1:nf
    nftri=nftri+1;
    tri(nftri,:)=[NN(1,i) NN(2,i) NN(3,i)];
    nftri=nftri+1;
    tri(nftri,:)=[NN(1,i) NN(3,i) NN(4,i)];
end
ligne=fgetl(fid);
ligne=fgetl(fid);
for i=1:nf
    ligne=fscanf(fid,'%g %g',6);
    xu(i)=ligne(1);
    yv(i)=ligne(2);
    zw(i)=ligne(3);
    u(i)=ligne(4);
    v(i)=ligne(5);
    w(i)=ligne(6);
end
status=fclose(fid);
figure;
trimesh(tri,x,y,z);
hold on;
quiver3(xu,yv,zw,u,v,w);
title('Mesh for Nemoh');
axis equal
savefig([projdir,filesep,'mesh',filesep,'axisym.fig'])
clear KH;
KH=zeros(6,6);
fid=fopen([projdir,filesep,'mesh',filesep,'KH.dat'],'r');
for i=1:6
    ligne=fscanf(fid,'%g %g',6);
    KH(i,:)=ligne;
end
status=fclose(fid);
clear XB YB ZB Mass WPA Inertia
Inertia=zeros(6,6);
fid=fopen([projdir,filesep,'mesh',filesep,'Hydrostatics.dat'],'r');
ligne=fscanf(fid,'%s',2);
XB=fscanf(fid,'%f',1);
ligne=fgetl(fid);
ligne=fscanf(fid,'%s',2);
YB=fscanf(fid,'%f',1);
ligne=fgetl(fid);
ligne=fscanf(fid,'%s',2);
ZB=fscanf(fid,'%f',1);
ligne=fgetl(fid);
ligne=fscanf(fid,'%s',2);
Mass=fscanf(fid,'%f',1)*1025.;
ligne=fgetl(fid);
ligne=fscanf(fid,'%s',2);
WPA=fscanf(fid,'%f',1);
status=fclose(fid);
clear ligne
fid=fopen([projdir,filesep,'mesh',filesep,'Inertia_hull.dat'],'r');
for i=1:3
    ligne=fscanf(fid,'%g %g',3);
    Inertia(i+3,4:6)=ligne;
end
Inertia(1,1)=Mass;
Inertia(2,2)=Mass;
Inertia(3,3)=Mass;

% Write Nemoh input file
fid=fopen([projdir,filesep,'Nemoh.cal'],'w');
fprintf(fid,'--- Environment ------------------------------------------------------------------------------------------------------------------ \n');
fprintf(fid,'1025.0				! RHO 			! KG/M**3 	! Fluid specific volume \n');
fprintf(fid,'9.81				! G			! M/S**2	! Gravity \n');
fprintf(fid,'%f                 ! DEPTH			! M		! Water depth\n',depth);
fprintf(fid,'0.	0.              ! XEFF YEFF		! M		! Wave measurement point\n');
fprintf(fid,'--- Description of floating bodies -----------------------------------------------------------------------------------------------\n');
fprintf(fid,'%g				! Number of bodies\n',nBodies);
for c=1:nBodies
    fprintf(fid,'--- Body %g -----------------------------------------------------------------------------------------------------------------------\n',c);
    %fprintf(fid,['axisym',int2str(c),'.dat\n']);
    fprintf(fid,['axisym.dat\n']);
    fprintf(fid,'%g %g			! Number of points and number of panels 	\n',nx(c),nf(c));
    fprintf(fid,'6				! Number of degrees of freedom\n');
    fprintf(fid,'1 1. 0.	0. 0. 0. 0.		! Surge\n');
    fprintf(fid,'1 0. 1.	0. 0. 0. 0.		! Sway\n');
    fprintf(fid,'1 0. 0. 1. 0. 0. 0.		! Heave\n');
    fprintf(fid,'2 1. 0. 0. %f %f %f		! Roll about a point\n',CoG(c,:));
    fprintf(fid,'2 0. 1. 0. %f %f %f		! Pitch about a point\n',CoG(c,:));
    fprintf(fid,'2 0. 0. 1. %f %f %f		! Yaw about a point\n',CoG(c,:));
    fprintf(fid,'6				! Number of resulting generalised forces\n');
    fprintf(fid,'1 1. 0.	0. 0. 0. 0.		! Force in x direction\n');
    fprintf(fid,'1 0. 1.	0. 0. 0. 0.		! Force in y direction\n');
    fprintf(fid,'1 0. 0. 1. 0. 0. 0.		! Force in z direction\n');
    fprintf(fid,'2 1. 0. 0. %f %f %f		! Moment force in x direction about a point\n',CoG(c,:));
    fprintf(fid,'2 0. 1. 0. %f %f %f		! Moment force in y direction about a point\n',CoG(c,:));
    fprintf(fid,'2 0. 0. 1. %f %f %f		! Moment force in z direction about a point\n',CoG(c,:));
    fprintf(fid,'0				! Number of lines of additional information \n');
end
fprintf(fid,'--- Load cases to be solved -------------------------------------------------------------------------------------------------------\n');
fprintf(fid,'%g %g	%f	%f	! Freq type 1,2,3=[rad/s,Hz,s],Number of wave frequencies, Min, and Max (rad/s)\n',1,length(wavefreq),wavefreq(1),wavefreq(end));
fprintf(fid,'%g	%f	%f		! Number of wave directions, Min and Max (degrees)\n',length(wavedir),wavedir(1),wavedir(end));
fprintf(fid,'--- Post processing ---------------------------------------------------------------------------------------------------------------\n');
fprintf(fid,'1	0.1	10.         ! IRF 				! IRF calculation (0 for no calculation), time step and duration\n');
fprintf(fid,'0                  ! Show pressure\n');
fprintf(fid,'0	0.	180.		! Kochin function 		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)\n');
fprintf(fid,'0	50	400. 400.   ! Free surface elevation 	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction\n');
fprintf(fid,'0                  ! Response Amplitude Operator (RAO), 0 no calculation, 1 calculated\n');
fprintf(fid,'1					! output freq type, 1,2,3=[rad/s,Hz,s]\n');
if QTFInput(1)==1
    fprintf(fid,'---QTF---\n');
    fprintf(fid,'%g         ! QTF flag, 1 is calculated \n',QTFInput(1));
    fprintf(fid,'%g	%f	%f	! Number of radial frequencies, Min, and Max values for the QTF computation \n',length(wavefreq),wavefreq(1),wavefreq(end));
    fprintf(fid,'%g         ! 0 Unidirection, Bidirection 1 \n', 0);
    fprintf(fid,'%g         ! Contrib, 1 DUOK, 2 DUOK+HASBO, 3 Full QTF (DUOK+HASBO+HASFS+ASYMP)\n',QTFInput(2));
    fprintf(fid,'NA 		! Name of free surface meshfile (Only for Contrib 3), type NA if not applicable \n');
    fprintf(fid,'0 	0	0	! Free surface QTF parameters: Re Nre NBessel (for Contrib 3)\n');
    fprintf(fid,'0          ! 1 Includes Hydrostatic terms of the quadratic first order motion, -[K]xi2_tilde \n');
    fprintf(fid,'1			! For QTFposProc, output freq type, 1,2,3=[rad/s,Hz,s]\n');
    fprintf(fid,'1         	! For QTFposProc, 1 includes DUOK in total QTFs, 0 otherwise\n');
    fprintf(fid,'1         	! For QTFposProc, 1 includes HASBO in total QTFs, 0 otherwise\n');
    fprintf(fid,'0         	! For QTFposProc, 1 includes HASFS+ASYMP in total QTFs, 0 otherwise\n');
else
    fprintf(fid,'---QTF---\n');
    fprintf(fid,'0         ! QTF flag, 1 is calculated \n');
end
fprintf(fid,'------\n');
status=fclose(fid);
fclose('all');

% if isfolder(projdir)
% disp('The same project name exists, ctrl-c for quit or enter for rewriting')
% pause;
% end
%system(['mv ','Mesh.cal ',projdir,filesep,'Mesh.cal']);

fid=fopen([projdir,filesep,'input_solver.txt'],'w');
fprintf(fid,'2				! Gauss quadrature (GQ) surface integration, N^2 GQ Nodes, specify N(1,4)\n');
fprintf(fid,'0.001			! eps_zmin for determine minimum z of flow and source points of panel, zmin=eps_zmin*body_diameter\n');
fprintf(fid,'1 				! 0 GAUSS ELIM.; 1 LU DECOMP.: 2 GMRES	!Linear system solver\n');
fprintf(fid,'10 1e-5 1000  	! Restart parameter, Relative Tolerance, max iter -> additional input for GMRES\n');
status=fclose(fid);

if QTFInput(1)==1
    system(['mkdir ',projdir,filesep,'Mechanics']);
    save([projdir,filesep,'Mechanics',filesep,'Inertia.dat'], 'Inertia','-ascii');
    save([projdir,filesep,'Mechanics',filesep,'Kh.dat'], 'KH','-ascii');
    Km=KH.*0;Badd=Km;
    save([projdir,filesep,'Mechanics',filesep,'Km.dat'], 'Km','-ascii');
    save([projdir,filesep,'Mechanics',filesep,'Badd.dat'], 'Badd','-ascii');
end

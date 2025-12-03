function [Idw,w,A,B,Fe,KH,Inertia,XB,YB,ZB,WPA] = BEM_Nemoh(projdir,ID_HydrosCal,ID_QTF,nBodies)

%% Nemoh HydroCalculation
if ID_QTF==1
    system(['mkdir ',projdir,filesep,'Motion']);
    system(['mkdir ',projdir,filesep,'results',filesep,'sources']);
end

% Calcul des coefficients hydrodynamiques
fprintf('\n------ Starting NEMOH ----------- \n');
system(['preProc ',projdir]);
fprintf('------ Solving BVPs ------------- \n');
system(['solver ',projdir]);
fprintf('------ Postprocessing results --- \n');
system(['postProc ',projdir]);



%% Lecture des resultats CA CM Fe
%clear Periode A B Famp Fphi Fe;

fid=fopen([projdir,filesep,'Nemoh.cal'],'r');
for i=1:6
    ligne=fgetl(fid);
end
nBodies=fscanf(fid,'%g',1);
for i=1:nBodies
    for ii=1:4
        ligne=fgetl(fid);
    end
    Ndof=fscanf(fid,'%g',1);
    for j=1:Ndof
        ligne=fgetl(fid);
    end
    ligne=fgetl(fid);
    Nforc=fscanf(fid,'%g',1);
    for j=1:Nforc
        ligne=fgetl(fid);
    end
    ligne=fgetl(fid);
end
ligne=fgetl(fid);
ligne=fgetl(fid);
ligne=fscanf(fid,'%g',2);
Idw=ligne(1);
nw=ligne(2);
fclose(fid);
fid=fopen([projdir,filesep,'results',filesep,'ExcitationForce.tec'],'r');
ligne=fgetl(fid);
for c=1:Nforc*nBodies
    ligne=fgetl(fid);
end
ligne=fgetl(fid);
for k=1:nw
    ligne=fscanf(fid,'%f',1+2*Nforc*nBodies);
    w(k)=ligne(1);
    for j=1:Nforc*nBodies
        Famp(k,j)=ligne(2*j);
        Fphi(k,j)=ligne(2*j+1);
    end
end
status=fclose(fid);
fid=fopen([projdir,filesep,'results',filesep,'RadiationCoefficients.tec'],'r');
ligne=fgetl(fid);
for i=1:Ndof*nBodies
    ligne=fgetl(fid);
end
for i=1:nBodies*Ndof
    ligne=fgetl(fid);
    for k=1:nw
        ligne=fscanf(fid,'%f',1+2*Ndof*nBodies);
        for j=1:Ndof*nBodies
            A(i,j,k)=ligne(2*j);
            B(i,j,k)=ligne(2*j+1);
        end
        ligne=fgetl(fid);
    end
end
status=fclose(fid);
% Expression des efforts d excitation de houle sous la forme Famp*exp(i*Fphi)
Fe=Famp.*(cos(Fphi)+1i*sin(Fphi));

end
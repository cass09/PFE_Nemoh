function [KH,Inertia,XB,YB,ZB,WPA] = hydrostatic_cal(projdir,nameNemohMesh,nemohInput)

isolatedBodyAnalysis = 1;
for c = 1:nemohInput.nBodies
    
    % Isolated body input nemoh datas
    isoNemohInput = nemohInput;
    isoNemohInput.nBodies = 1;
    isoNemohInput.CG = nemohInput.CG(c,:);
    isoNemohInput.transMatrix = nemohInput.transMatrix(c,:);
    isoNemohInput.nodes = nemohInput.nodes(c,:);
    isoNemohInput.panels = nemohInput.panels(c,:);
    
    %% Writing input files
    writeNemohCalFile(projdir,nameNemohMesh,isoNemohInput);   % Write Nemoh.cal
    tX = isoNemohInput.transMatrix(1,1); tY = isoNemohInput.transMatrix(1,2);
    % writeMeshCalFile(projdir,nameNemohMesh(c),isolatedBodyAnalysis,tX,tY,isoNemohInput.CG,nemohInput.panels(c));  % Write Nemoh.cal
    
    
    system(['mkdir ',projdir,filesep,'Mechanics']);
    %% Calculation
    fprintf('\n------ Starting NEMOH ----------- \n');
    system(['preProc ',projdir]);
    fprintf('------ computes Hydrostatic ------------- \n');
    system(['hydrosCal ',projdir]);


%% Reading hydrostatic properties from files
    % Hydrostatic matrix KH
    fid=fopen([projdir,filesep,'mesh',filesep,'KH.dat'],'r');
    for i=1:6
        ligne=fscanf(fid,'%g %g',6);
        KH(c,i,:)=ligne;
    end
    status=fclose(fid);

    % CoB, Mass, WPA
    fid=fopen([projdir,filesep,'mesh',filesep,'Hydrostatics.dat'],'r');
    ligne=fscanf(fid,'%s',2);
    XB(c)=fscanf(fid,'%f',1);
    ligne=fgetl(fid);
    ligne=fscanf(fid,'%s',2);
    YB(c)=fscanf(fid,'%f',1);
    ligne=fgetl(fid);
    ligne=fscanf(fid,'%s',2);
    % ZB(c)=fscanf(fid,'%f',1);
    str = fscanf(fid,'%s',1);
    val = str2double(str);
    if isempty(str) || all(str=='*') || isnan(val)
        val = 0;
    end
    ZB(c) = val;
    ligne=fgetl(fid);
    ligne=fscanf(fid,'%s',2);
    Mass(c)=fscanf(fid,'%f',1)*1025.;
    ligne=fgetl(fid);
    ligne=fscanf(fid,'%s',3);
    WPA(c)=fscanf(fid,'%f',1);
    status=fclose(fid);
    clear ligne
    % Inertia
    fid=fopen([projdir,filesep,'mesh',filesep,'Inertia_hull.dat'],'r');
    for i=1:3
        ligne=fscanf(fid,'%g %g',3);
        Inertia(c,i+3,4:6)=ligne;
    end
    Inertia(c,1,1)=Mass(c);
    Inertia(c,2,2)=Mass(c);
    Inertia(c,3,3)=Mass(c);


end

end
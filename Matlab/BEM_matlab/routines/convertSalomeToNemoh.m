function convertSalomeToNemoh(folderSalome,folderNemoh,nameSalomeMesh,nameNemohMesh)
%% Conversion of the mesh file format from SALOME to Nemoh 
%salome-nemoh conversion format


    result = dlmread([folderSalome,filesep,nameSalomeMesh,'.dat']);
    fid = fopen ([folderNemoh,filesep,nameNemohMesh,'.dat'], 'w');
    fprintf(fid, '%g %g\n',[2 0]);
    
    % Translation x coord
    result(2:result(1,1)+1,2) = result(2:result(1,1)+1,2);
    % Translation y coord
    result(2:result(1,1)+1,3) = result(2:result(1,1)+1,3);
    % Translation z coord
    result(2:result(1,1)+1,4) = result(2:result(1,1)+1,4);

    for jj = 1:result(1,1)
       fprintf(fid,'%g %f %f %f\n',result(jj+1,1:4));   % write nodes coordinates  
    end
    
    % Panels 
    rows102 = find(result(:,2)==102); %index where the connectivity table refers to 2 nodes edge reference 
    rows203 = find(result(:,2)==203); %index where the connectivity table refers to 3 nodes triangle reference
    rows204 = find(result(:,2)==204); %index where the connectivity table refers to 4 nodes quadrangular reference
    
    fprintf(fid, '%g %g %g %g\n',[0 0 0 0]);  % delimeter row
    
    for kk = 1:length(rows203)
       fprintf(fid,'%g %g %g %g\n',result(rows203(kk),3),result(rows203(kk),4),result(rows203(kk),5),result(rows203(kk),3));    
    end
    
    for kk = 1:length(rows204)
       fprintf(fid,'%g %g %g %g\n',result(rows204(kk),3),result(rows204(kk),4),result(rows204(kk),5),result(rows204(kk),6));     
    end
    
    fprintf(fid, '%g %g %g %g\n',[0 0 0 0]);
    
    fclose(fid)

    nodes = result(1,1);
    panels= length(rows203)+length(rows204);
    
end
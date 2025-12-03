function [nodes,panels]= translateMesh(nBodies,folderNemoh,nameNemohMesh,nameNemohMesh_trasl,transVector)
%% Conversion of the mesh file format from SALOME to Nemoh 
%salome-nemoh conversion format


    result = dlmread([folderNemoh,filesep,nameNemohMesh,'.dat']);
    fid = fopen ([folderNemoh,filesep,nameNemohMesh_trasl,'.dat'], 'w');
    fprintf(fid, '%g %g\n',[2 result(1,2)]);
    %
    row_all_zeros = find(all(result == 0,2));
    % number of nodes
    nodes = row_all_zeros(1) -2; % not count first row and zeros row delimeter
    panels = row_all_zeros(2) -row_all_zeros(1) - 1; %  not count last zeros row delimeter
   
    % Translation x coord
    result(2:nodes+1,2) = result(2:nodes+1,2)+transVector(1);
    % Translation y coord
    result(2:nodes+1,3) = result(2:nodes+1,3)+transVector(2);
    % Translation z coord
    result(2:nodes+1,4) = result(2:nodes+1,4)+transVector(3);
   
    % print Nodes
    for jj = 1:nodes
       fprintf(fid,'%g %f %f %f\n',result(jj+1,1:4));   % write nodes coordinates  
    end
    
    % print Panels 
    fprintf(fid, '%g %g %g %g\n',[0 0 0 0]);  % delimeter row

    for kk = 1:row_all_zeros(2)-row_all_zeros(1)
       fprintf(fid,'%g %g %g %g\n',result(row_all_zeros(1)+kk,:));    
    end

    fclose(fid)

    %% DISPLAY MESH
% Nodes coordinates
x = result(2:row_all_zeros(1)-1,2);
y = result(2:row_all_zeros(1)-1,3);
z = result(2:row_all_zeros(1)-1,4);

% Panels 
NN(1,:) = result(row_all_zeros(1)+1:end-1,1);
NN(2,:) = result(row_all_zeros(1)+1:end-1,2);
NN(3,:) = result(row_all_zeros(1)+1:end-1,3);
NN(4,:) = result(row_all_zeros(1)+1:end-1,4);

    nftri=0;
    for i=1:panels
        nftri=nftri+1;
        tri(nftri,:)=[NN(1,i) NN(2,i) NN(3,i)];
        nftri=nftri+1;
        tri(nftri,:)=[NN(1,i) NN(3,i) NN(4,i)];
    end



% figure(199);
% hold on
% trimesh(tri,x,y,z);
% axis equal
% xlabel('x (m)'), ylabel('y (m)'), zlabel('z (m)')
% view([-45 25])



end
function writeMeshCalFile(projdir,nameNemohMesh,nBodies,tX,tY,CG,nfobj)
%   Creation des fichiers de calcul du maillage


for c=1:nBodies
    fprintf('\n -> Meshing body number %g \n',c);

    fid=fopen([projdir,filesep,'Mesh.cal'],'w');  %
    fprintf(fid,[nameNemohMesh{c,1},'\n'],1);
    fprintf(fid,'0 \n %f %f \n ',tX(c), tY(c));   % 1 if a half symmetric body mesh, about (xOz)
    fprintf(fid,'%f %f %f \n',CG(c,:));    
    fprintf(fid,'%g \n 2 \n 0. \n 1.\n',nfobj(c));
    fprintf(fid,'%f \n',1025.);
    fprintf(fid,'%f \n',9.81);
    status=fclose(fid);
end
end


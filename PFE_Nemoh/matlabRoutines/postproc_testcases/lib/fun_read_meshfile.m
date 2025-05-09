function [vert,conn]=fun_read_meshfile(meshFile,NNodes,NPanels)
fid=fopen(meshFile,'r');
ligne=fgetl(fid);
vert=zeros(NNodes,3);
conn=zeros(NPanels,4);
for inode=1:NNodes
    val=fscanf(fid,'%g',4);
    vert(inode,1:3)=val(2:4);
    ligne=fgetl(fid);
end
ligne=fgetl(fid);
for ipanel=1:NPanels
    val=fscanf(fid,'%d',4);
    conn(ipanel,1:4)=val(1:4);
    ligne=fgetl(fid);
end
fclose(fid);
end
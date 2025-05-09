function dat=fun_readHydrostarOutput(filename,nline0)

fid=fopen(filename,'r');
for ii=1:nline0
    ligne=fgetl(fid);
end
iter=1;
IdFlag=1;
while IdFlag==1
    val=fscanf(fid,'%g',3);
    if ~isempty(val)
        dat(iter,1:3)=val;
        ligne=fgetl(fid);
        iter=iter+1;
    else
        IdFlag=0;
    end
    
end
fclose(fid);
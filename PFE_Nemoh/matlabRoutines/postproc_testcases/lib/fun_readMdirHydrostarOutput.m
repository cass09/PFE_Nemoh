function dat=fun_readMdirHydrostarOutput(filename,nline0,Idbeta,Nbetadat)

fid=fopen(filename,'r');
for ii=1:nline0
    ligne=fgetl(fid);
end
iter=1;
IdFlag=1;
while IdFlag==1
    val=fscanf(fid,'%g',Nbetadat*2+1);
    if ~isempty(val)
        dat(iter,1)=val(1);
        dat(iter,2)=val(Idbeta+1);
        dat(iter,3)=val(Nbetadat+Idbeta+1);
        ligne=fgetl(fid);
        iter=iter+1;
    else
        IdFlag=0;
    end
    
end
fclose(fid);
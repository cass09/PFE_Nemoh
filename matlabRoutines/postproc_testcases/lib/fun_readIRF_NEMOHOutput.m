function [IRFdat]=fun_readIRF_NEMOHOutput(rep)

    fid=fopen([rep,filesep,'Nemoh.cal'],'r');
    for i=1:6
        ligne=fgetl(fid);
    end
    nBodies=fscanf(fid,'%g',1);
    ligne=fgetl(fid);
    for i=1:nBodies
        for j=1:3
            ligne=fgetl(fid);
        end
        Ndof=fscanf(fid,'%g',1);
        ligne=fgetl(fid);
        for j=1:Ndof
            ligne=fgetl(fid);
        end
        Nforce=fscanf(fid,'%g',1);
        ligne=fgetl(fid);
        for j=1:Nforce
            ligne=fgetl(fid);
        end
        Nadd=fscanf(fid,'%g',1);
        ligne=fgetl(fid);
        if Nadd>0
            Ndeflm=fscanf(fid,'%g',1);
            ligne=fgetl(fid);
            for j=1:Nadd-1
                ligne=fgetl(fid);
            end
        else
            Ndeflm=0;
        end
    end
    ligne=fgetl(fid);
    ligne=fscanf(fid,'%d',2);
    nw=ligne(2);
    ligne=fgetl(fid);
    ligne=fscanf(fid,'%g',3);
    Nbeta=ligne(1);
    betamin=ligne(2);
    betamax=ligne(3);
    ligne=fgetl(fid);
    ligne=fgetl(fid);
    ligne=fscanf(fid,'%d %f %f',3);
    IRFswitch=ligne(1);
    IRFNt=2*(ligne(3)/ligne(2))-1;
    fclose(fid);
    beta=linspace(betamin,betamax,Nbeta);
    
    IRFdat=zeros(IRFNt,1+Nforce*nBodies,Nbeta);
    ww=zeros(nw,1);
   
    fid=fopen([rep,filesep,'results',filesep,'IRF_excForce.dat'],'r');
    ligne=fgetl(fid);
    Nforce=Nforce+Ndeflm;
    Ndof=Ndof+Ndeflm;
   
    for ibeta=1:Nbeta
        ligne=fgetl(fid);
        for k=1:IRFNt
                IRFdat(k,:,ibeta)=fscanf(fid,'%f',1+Nforce*nBodies);
                ligne=fgetl(fid);
        end
    end
    status=fclose(fid);

    
   
end
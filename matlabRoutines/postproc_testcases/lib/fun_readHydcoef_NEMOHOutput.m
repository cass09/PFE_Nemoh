function [ww,Fdat,AmatNem,BmatNem,RAONem]=fun_readHydcoef_NEMOHOutput(rep,IDF,IDABMAT,IDRAO)

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
    fclose(fid);
    beta=linspace(betamin,betamax,Nbeta);

    Fdat=zeros(nw,2*Nforce*nBodies,Nbeta);
    ww=zeros(nw,1);
    if (IDF==1)
    fid=fopen([rep,filesep,'results',filesep,'ExcitationForce.tec'],'r');
    ligne=fgetl(fid);
    Nforce=Nforce+Ndeflm;
    Ndof=Ndof+Ndeflm;
    for c=1:Nforce*nBodies
        ligne=fgetl(fid);
    end
    for ibeta=1:Nbeta
        ligne=fgetl(fid);
        for k=1:nw
                ligne=fscanf(fid,'%f',1+2*Nforce*nBodies);
                ww(k)=ligne(1);
                for IDDOF=1:Nforce*nBodies
                Fdat(k,IDDOF,ibeta)=ligne(2*IDDOF);
                Fdat(k,Nforce*nBodies+IDDOF,ibeta)=ligne(2*IDDOF+1);
                end
               ligne=fgetl(fid);
        end
    end
    status=fclose(fid);
    end

    AmatNem=zeros(nw,nBodies*Ndof,nBodies*Ndof);
    BmatNem=zeros(nw,nBodies*Ndof,nBodies*Ndof);
    if (IDABMAT==1)
    fid=fopen([rep,filesep,'results',filesep,'RadiationCoefficients.tec'],'r');
    ligne=fgetl(fid);
    for i=1:Ndof*nBodies
        ligne=fgetl(fid);
    end

    for i=1:nBodies*Ndof
        ligne=fgetl(fid);
        for k=1:nw
            ligne=fscanf(fid,'%f',1+2*Ndof*nBodies);
            ww(k)=ligne(1);
            for j=1:i
                AmatNem(k,i,j)=ligne(2*j);
                BmatNem(k,i,j)=ligne(2*j+1);
            end
            ligne=fgetl(fid);
        end
    end
    status=fclose(fid);
    end

%     % Expression des efforts d excitation de houle sous la forme Famp*exp(i*Fphi)
%     Fe=Famp.*(cos(Fphi)+1i*sin(Fphi));
RAONem=zeros(nw,2*Ndof*nBodies,Nbeta);
if IDRAO==1
        fid=fopen([rep,filesep,'Motion',filesep,'RAO.dat'],'r');

        ligne=fgetl(fid);
        for ibeta=1:Nbeta
            ligne=fgetl(fid);
            for iw=1:nw
                ligne=fscanf(fid,'%g',1+2*Ndof*nBodies);
                ww(iw)=ligne(1);
                for IDDOF=1:Ndof*nBodies
                RAONem(iw,IDDOF,ibeta)=ligne(1+IDDOF);
                RAONem(iw,Ndof*nBodies+IDDOF,ibeta)=deg2rad(ligne(IDDOF+1+Ndof*nBodies));
                end
                ligne=fgetl(fid);
            end
        end
        status=fclose(fid);
end
end
function [QTFTotRe,QTFTotIm,QTFTotMod,ww]=fun_readNEMOHQTFOutput_N...
        (filename,nw,dof_data,DOF,betaId,Nbetadat,qtftype,SwitchBiDir)
IdhalfdiagData=1; %depend on the NEMOH data

fid=fopen(filename,'r');
nline0=1;
for ii=1:nline0
    ligne=fgetl(fid);
end

ww=zeros(1,nw);
QTFTotRe=zeros(nw,nw);
QTFTotIm=zeros(nw,nw);
QTFTotMod=zeros(nw,nw);
Nbeta2dat=Nbetadat;
if (SwitchBiDir==0),Nbeta2dat=1;end
for IDdof=1:6
    if (dof_data(IDdof)==0), continue; end
    for ib1=1:Nbetadat
        for ib2=1:Nbeta2dat   
            for ii=1:nw
                for jj=1:ii
                    val=fscanf(fid,'%f %f %f %f %d %g %g %g %g',9);
                    if (round(val(3))==betaId(1)&&round(val(4))==betaId(2))
                        if (IDdof==val(5)&& IDdof==DOF)
                            if (ii==1)
                            disp(['DOF=',num2str(IDdof)])
                            disp(['NEMOH: beta1= ',num2str(round(val(3))),...
                                ', beta2= ',num2str(round(val(4)))])
                            end
                            QTFTotRe(ii,jj)=val(8);
                            QTFTotIm(ii,jj)=val(9);
                            QTFTotMod(ii,jj)=val(6);
                        end
                    end
                end
                ww(ii)=val(1);
            end
        end
    end
end
fclose(fid);

QTFTotRe=tril(QTFTotRe)+triu(QTFTotRe.',1);
QTFTotMod=tril(QTFTotMod)+triu(QTFTotMod.',1);

if strcmp(qtftype,'M')
QTFTotIm=(tril(QTFTotIm)-triu(QTFTotIm.',1));
else
QTFTotIm=(tril(QTFTotIm)+triu(QTFTotIm.',1));
end
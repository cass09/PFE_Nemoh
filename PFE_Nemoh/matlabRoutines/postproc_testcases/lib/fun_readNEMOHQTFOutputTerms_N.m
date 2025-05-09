function [QTFTotRe,QTFTotIm,QTFTotMod,ww]=fun_readNEMOHQTFOutputTerms_N...
    (filename,nw,dof_data,DOF,betaId,Nbetadat,qtftype,rho,g)

fid=fopen(filename,'r');
nline0=1;
for ii=1:nline0
    ligne=fgetl(fid);
end

ww=zeros(1,nw);
QTFTotRe=zeros(nw,nw);
QTFTotIm=zeros(nw,nw);
QTFTotMod=zeros(nw,nw);
for ib1=1:Nbetadat
    for ib2=1:Nbetadat
        for i0=0:nw-1
            for ii=i0+1:nw
                jj=ii-i0;
                for IDdof=1:6
                    if (dof_data(IDdof)==0), continue; end
                    val=fscanf(fid,'%f %f %f %f %d %g %g',7);
                    if (round(rad2deg(val(3)))==betaId(1)&&round(rad2deg(val(4)))==betaId(2))
                        if (IDdof==val(5)&& IDdof==DOF)
                            if (ii==1 && jj==1)
                                disp(['DOF=',num2str(IDdof)])
                                disp(['NEMOH: beta1= ',num2str(round(rad2deg(val(3)))),...
                                    ', beta2= ',num2str(round(rad2deg(val(4))))])
                            end
                            QTFTotRe(ii,jj)=val(6)/rho/g;
                            QTFTotIm(ii,jj)=val(7)/rho/g;
                            QTFTotMod(ii,jj)=sqrt(val(6)^2+val(7)^2)/rho/g;
                            
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
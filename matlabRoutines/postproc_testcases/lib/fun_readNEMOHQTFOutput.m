function [QTFTotRe,QTFTotIm,QTFTotMod,ww]=fun_readNEMOHQTFOutput(filename,nw,dof_data,DOF,qtftype)
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
for IDdof=1:6
    if (dof_data(IDdof)==0), continue; end
    for ii=1:nw
        for jj=ii^IdhalfdiagData:nw
            val=fscanf(fid,'%f %f %f %f %d %g %g %g %g',9);
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
        ww(ii)=val(1);
    end
end
fclose(fid);

QTFTotRe=triu(QTFTotRe)+tril(QTFTotRe.',-1);
QTFTotMod=triu(QTFTotMod)+tril(QTFTotMod.',-1);

if strcmp(qtftype,'M')
QTFTotIm=(triu(QTFTotIm)-tril(QTFTotIm.',-1));
else
QTFTotIm=(triu(QTFTotIm)+tril(QTFTotIm.',-1));
end
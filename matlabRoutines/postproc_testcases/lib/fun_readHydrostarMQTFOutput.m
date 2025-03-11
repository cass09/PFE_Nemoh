function [QTFTotRe,QTFTotIm,ww]=fun_readHydrostarMQTFOutput(filename,nline0,nw,betaId,nbeta,qtftype)

fid=fopen(filename,'r');
for ii=1:nline0-1
    ligne=fgetl(fid);
end
iter=1;
IdFlag=1;
datRe=zeros(nw,nw+1);
datIm=zeros(nw,nw+1);
ww=zeros(nw,1);

for ib1=1:nbeta
    val= fscanf(fid, '%s',3);
    beta1=str2num(val(12:15));
    for ib2=1:nbeta
        val= fscanf(fid, '%s',3);
        beta2=str2num(val(11:14));
        if (beta1==betaId(1) && beta2==betaId(2))
             disp(['HYDROSTAR: beta1= ',num2str(beta1),...
                               ', beta2= ',num2str(beta2)])
        end
        for ii=1:nw
            val=fscanf(fid,'%g',2*(nw+1)+1);
            if ~isempty(val)
                if (beta1==betaId(1) && beta2==betaId(2))
                    datRe(ii,1:end)=val(2:2:end-1);
                    datIm(ii,1:end)=val(3:2:end);
                    ww(ii)=val(1);
                    ligne=fgetl(fid);
                end
            else
                IdFlag=0;
            end
        end
        ligne=fgetl(fid);
        ligne=fgetl(fid);
    end
end


% assignin('base','datIm',datIm);

dw=zeros(1,nw+1);
dw(2:end)=ww.';

fclose(fid);
QTFTotRe=zeros(nw,nw);
QTFTotIm=zeros(nw,nw);

for ii=1:nw
    for jj=ii:nw
        inddw=closest(dw,ww(jj)-ww(ii));
        QTFTotRe(ii,jj)=datRe(ii,inddw);
        QTFTotIm(ii,jj)=datIm(ii,inddw);
    end
end
QTFTotRe=triu(QTFTotRe)+tril(QTFTotRe.',-1);
%due to different potential converence with NEMOH the imaginary part of
%Hydrostar adjustes
%the diagonal part is usually 0 for the imaginar except for bidirectional
%case.
if strcmp(qtftype,'M')
QTFTotIm=QTFTotIm.*eye(size(QTFTotIm))-(triu(QTFTotIm,1)-tril(QTFTotIm.',-1));
else
QTFTotIm=-(triu(QTFTotIm)+tril(QTFTotIm.',-1));
end
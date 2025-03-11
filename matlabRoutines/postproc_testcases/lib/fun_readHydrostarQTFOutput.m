function [QTFTotRe,QTFTotIm,ww]=fun_readHydrostarQTFOutput(filename,nline0,nw,qtftype)

fid=fopen(filename,'r');
for ii=1:nline0
    ligne=fgetl(fid);
end
iter=1;
IdFlag=1;
datRe=zeros(nw,nw+1);
datIm=zeros(nw,nw+1);
ww=zeros(nw,1);
for ii=1:nw
    val=fscanf(fid,'%g',2*(nw+1)+1);
    if ~isempty(val)
        datRe(ii,1:end)=val(2:2:end-1);
        datIm(ii,1:end)=val(3:2:end); 
        ww(ii)=val(1);
        ligne=fgetl(fid);
    else
        IdFlag=0;
    end   
end
% assignin('base','datRe',datRe);
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
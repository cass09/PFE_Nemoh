function funP_contour(axesfig,XX,YY,ZZ,N,Zlev,col,lw)
tcon = hgtransform('Parent',axesfig);
if ~isempty(N)
[c,h]=contour(XX,YY,ZZ,N,'parent',tcon);
else
[c,h]=contour(XX,YY,ZZ,'parent',tcon);    
end
clabel(c,h,'FontSize',10,'FontWeight','bold','Color','k')
if ~isempty(col)
set(h,'linecolor',col,'linewidth',lw)    
end
set(tcon, 'Matrix',makehgtform('translate',[0 0 Zlev]));
%plot properties     
fontsize=15;
linewidth=2;
lineobj = findobj('type', 'line');
set(lineobj, 'linewidth', linewidth);
%set(lineobj, 'linestyle', '--');
textobj = findobj('type', 'text'); %buat legend
%set(textobj, 'fontunits', 'points');
set(textobj, 'fontsize', fontsize);
set(textobj, 'fontweight', 'bold');
xl = get(gca, 'xlabel');
set(xl,'fontweight', 'bold');
set(xl,'fontsize', fontsize);
yl = get(gca, 'ylabel');
set(yl,'fontweight', 'bold');
set(yl,'fontsize', fontsize);
set(gca, 'FontSize', fontsize);
set(gca, 'Fontweight', 'bold');
tit = get(gca, 'title');
set(tit, 'FontSize', fontsize);
set(tit, 'Fontweight', 'bold');
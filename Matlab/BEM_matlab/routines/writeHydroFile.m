function writeHydroFile(projdir,nBodies,CG)
fid=fopen([projdir,filesep,'Hydro.txt'],'w');  %
fprintf(fid,'%g ! Number of bodies \n',nBodies);
for c=1:nBodies
    fprintf(fid,'%f %f %f \n',CG(c,:));    
status=fclose(fid);

end
end


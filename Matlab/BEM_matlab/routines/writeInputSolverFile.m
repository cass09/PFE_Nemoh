function writeInputSolverFile(projdir)
fid=fopen([projdir,filesep,'input_solver.txt'],'w');
fprintf(fid,'2				! Gauss quadrature (GQ) surface integration, N^2 GQ Nodes, specify N(1,4)\n');
fprintf(fid,'0.001			! eps_zmin for determine minimum z of flow and source points of panel, zmin=eps_zmin*body_diameter\n');
fprintf(fid,'0 				! 0 GAUSS ELIM.; 1 LU DECOMP.: 2 GMRES	!Linear system solver\n');
fprintf(fid,'10 1e-5 1000  	! Restart parameter, Relative Tolerance, max iter -> additional input for GMRES\n');
status=fclose(fid);
end
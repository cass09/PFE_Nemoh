function writeNemohCalFile(projdir,namemeshNemoh,nemohInput)
% Write Nemoh input file

fid=fopen([projdir,filesep,'Nemoh.cal'],'w');
fprintf(fid,'--- Environment ------------------------------------------------------------------------------------------------------------------ \n');
fprintf(fid,'%g				! RHO 			! KG/M**3 	! Fluid specific volume \n',nemohInput.rho);
fprintf(fid,'%f				! G			! M/S**2	! Gravity \n',nemohInput.g);
fprintf(fid,'%f                 ! DEPTH			! M		! Water depth\n',nemohInput.depth);
fprintf(fid,'0.	0.              ! XEFF YEFF		! M		! Wave measurement point\n');
fprintf(fid,'--- Description of floating bodies -----------------------------------------------------------------------------------------------\n');
fprintf(fid,'%g				! Number of bodies\n',nemohInput.nBodies);
for c=1:nemohInput.nBodies
    fprintf(fid,'--- Body %g -----------------------------------------------------------------------------------------------------------------------\n',c);
    fprintf(fid,[namemeshNemoh{c,1},'.dat\n']);
    fprintf(fid,'%g %g			! Number of points and number of panels 	\n',nemohInput.nodes(c),nemohInput.panels(c));
    fprintf(fid,'%g				! Number of degrees of freedom\n',nemohInput.nDegFreedom);
    if nemohInput.DOF(1)==1
        fprintf(fid,'1 1. 0.	0. 0. 0. 0.		! Surge\n');
    end
    if nemohInput.DOF(2)==1
        fprintf(fid,'1 0. 1.	0. 0. 0. 0.		! Sway\n');
    end
    if nemohInput.DOF(3)==1
        fprintf(fid,'1 0. 0. 1. 0. 0. 0.		! Heave\n');
    end
    if nemohInput.DOF(4)==1
        fprintf(fid,'2 1. 0. 0. %f %f %f		! Roll about a point\n',nemohInput.CG(c,:));
    end
    if nemohInput.DOF(5)==1
        fprintf(fid,'2 0. 1. 0. %f %f %f		! Pitch about a point\n',nemohInput.CG(c,:));
    end
    if nemohInput.DOF(6)==1
        fprintf(fid,'2 0. 0. 1. %f %f %f		! Yaw about a point\n',nemohInput.CG(c,:));
    end
    fprintf(fid,'%g				! Number of resulting generalised forces\n',nemohInput.nDegFreedom);
    if nemohInput.DOF(1)==1
        fprintf(fid,'1 1. 0.	0. 0. 0. 0.		! Force in x direction\n');
    end
    if nemohInput.DOF(2)==1
        fprintf(fid,'1 0. 1.	0. 0. 0. 0.		! Force in y direction\n');
    end
    if nemohInput.DOF(3)==1
        fprintf(fid,'1 0. 0. 1. 0. 0. 0.		! Force in z direction\n');
    end
    if nemohInput.DOF(4)==1
        fprintf(fid,'2 1. 0. 0. %f %f %f		! Moment force in x direction about a point\n',nemohInput.CG(c,:));
    end
    if nemohInput.DOF(5)==1
        fprintf(fid,'2 0. 1. 0. %f %f %f		! Moment force in y direction about a point\n',nemohInput.CG(c,:));
    end
    if nemohInput.DOF(6)==1
        fprintf(fid,'2 0. 0. 1. %f %f %f		! Moment force in z direction about a point\n',nemohInput.CG(c,:));
    end
    fprintf(fid,'0				! Number of lines of additional information \n');
end
fprintf(fid,'--- Load cases to be solved -------------------------------------------------------------------------------------------------------\n');
fprintf(fid,'%g %g	%f	%f	! Freq type 1,2,3=[rad/s,Hz,s],Number of wave frequencies, Min, and Max (rad/s)\n',nemohInput.wavetype,length(nemohInput.wavefreq),nemohInput.wavefreq(1),nemohInput.wavefreq(end));
fprintf(fid,'%g	%f	%f		! Number of wave directions, Min and Max (degrees)\n',nemohInput.wavedir(1),nemohInput.wavedir(2),nemohInput.wavedir(3));
fprintf(fid,'--- Post processing ---------------------------------------------------------------------------------------------------------------\n');
fprintf(fid,'%g	%f %f         ! IRF 				! IRF calculation (0 for no calculation), time step and duration\n',nemohInput.IRFcalculation(1),nemohInput.IRFcalculation(2),nemohInput.IRFcalculation(3));
fprintf(fid,'%g                  ! Show pressure\n',nemohInput.showpressure);
fprintf(fid,'%g	%f %f		! Kochin function 		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)\n',nemohInput.kochinfunction(1),nemohInput.kochinfunction(2),nemohInput.kochinfunction(3));
fprintf(fid,'%g	%g %f %f   ! Free surface elevation 	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction\n',nemohInput.freesurfaceelevation(1),nemohInput.freesurfaceelevation(2),nemohInput.freesurfaceelevation(3),nemohInput.freesurfaceelevation(4));
fprintf(fid,'%g                  ! Response Amplitude Operator (RAO), 0 no calculation, 1 calculated\n', nemohInput.RAO);
fprintf(fid,'1					! output freq type, 1,2,3=[rad/s,Hz,s]\n');
fprintf(fid,'--- Interaction Theory ---------------------------------------------------------------------------------------------------------------\n');
fprintf(fid,'%f %g %g		! Interaction Theory Cylindrical Envelop\n', nemohInput.OuterCylinder(1), nemohInput.OuterCylinder(2), nemohInput.OuterCylinder(3));
fprintf(fid,'%g 0 		! Interaction Theory with sources\n', nemohInput.Source);
fprintf(fid,'0			! run IT\n');
if nemohInput.QTFInput(1)==1
    fprintf(fid,'---QTF---\n');
    fprintf(fid,'%g         ! QTF flag, 1 is calculated \n',nemohInput.QTFInput(1));
    fprintf(fid,'%g	%f	%f	! Number of radial frequencies, Min, and Max values for the QTF computation \n',length(nemohInput.wavefreq),nemohInput.wavefreq(1),nemohInput.wavefreq(end));
    fprintf(fid,'%g         ! 0 Unidirection, Bidirection 1 \n', 0);
    fprintf(fid,'%g         ! Contrib, 1 DUOK, 2 DUOK+HASBO, 3 Full QTF (DUOK+HASBO+HASFS+ASYMP)\n',QTFInput(2));
    fprintf(fid,'NA 		! Name of free surface meshfile (Only for Contrib 3), type NA if not applicable \n');
    fprintf(fid,'0 	0	0	! Free surface QTF parameters: Re Nre NBessel (for Contrib 3)\n');
    fprintf(fid,'0          ! 1 Includes Hydrostatic terms of the quadratic first order motion, -[K]xi2_tilde \n');
    fprintf(fid,'1			! For QTFposProc, output freq type, 1,2,3=[rad/s,Hz,s]\n');
    fprintf(fid,'1         	! For QTFposProc, 1 includes DUOK in total QTFs, 0 otherwise\n');
    fprintf(fid,'1         	! For QTFposProc, 1 includes HASBO in total QTFs, 0 otherwise\n');
    fprintf(fid,'0         	! For QTFposProc, 1 includes HASFS+ASYMP in total QTFs, 0 otherwise\n');
else
    fprintf(fid,'---QTF---\n');
    fprintf(fid,'0         ! QTF flag, 1 is calculated \n');
end
fprintf(fid,'------\n');
status=fclose(fid);
fclose('all');
end
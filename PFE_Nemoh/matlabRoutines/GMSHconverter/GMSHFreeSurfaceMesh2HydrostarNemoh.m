%FS Mesh converter from GMSH to NEMOH and HYDROSTAR
%Created by Ruddy Kurnia, ECN-LHEAA 2021-2022
%Note only quadrilateral mesh allow

clear all;clc;
OutputFile=1;%1 for Nemoh, 2 for Hydrostar
IncludeBdyLine=1;    %1 included, 0 not.. Only for NEMOH
Meshfile='OC5CONFIG_MAIN_FREESURFACE_R50'; % input file name without ext
MeshfileOutname='OC5CONFIG_MAIN_FREESURFACE_R50'; % output file name without ext
dir='C:\Users\rkurnia\Documents\Test-VirtualBox\gmsh-4.8.4-Linux64\bin\OC5_CONFIG_MAIN\FS_mesh\';

% Meshfile='OC5CONFIG1_FREESURFACE_R550'; % input file name without ext
% MeshfileOutname='OC5CONFIG1_FREESURFACE_R550'; % output file name without ext
% dir='.\OC5_CONFIG_meshes\freesurface_meshes';

fid=fopen([dir,filesep,Meshfile,'.msh'],'r');
if OutputFile==1
IdSym=input('\n -input 1 if half-body mesh, 0 otherwise: ');
fidWrite=fopen([dir,filesep,MeshfileOutname,'_FSNem.dat'],'w');
fprintf(fidWrite,'%i %s %s %s %s %s %s %s\n',IdSym,'[NPoints]', '[NPanels]', '[NBdyLines]', 'to be specified'); %%  0 is not symmetric, 1 is symmetric  
elseif OutputFile==2
disp('Number of point of the body panel is needed as a start point ID of the free surface point ')
NPstarts=input('\n -input number of point of Body panel: ');
fidWrite=fopen([dir,filesep,MeshfileOutname,'_FSHS.dat'],'w');
fprintf(fidWrite,'PROJECT xxxxxxxxxxxxxx \n');
fprintf(fidWrite,'USER    xxxxxxxxxxxxxx \n');
fprintf(fidWrite,'\n');
end
ligne='****';
while ~strcmpi(ligne,'$Nodes') 
ligne=fgetl(fid);
end
val=fscanf(fid,'%g',4);
NNodes=val(2);
XYZ=zeros(NNodes,3);
ligne=fgetl(fid);
iter=1;



while iter<=NNodes
 val=fscanf(fid,'%g',4); 
 ligne=fgetl(fid);  
 Nnodetemp=val(4);  
 for kk=1:Nnodetemp
 ligne=fgetl(fid);      
 end
 for kk=1:Nnodetemp
 val=fscanf(fid,'%g',3);
 XYZ(iter,:)=val.';
 if OutputFile==1
 fprintf(fidWrite,'%i %g %g	%g\n',iter,val(1),val(2),val(3));   
 else
 fprintf(fidWrite,'%i %g %g	%g\n',iter+NPstarts,val(1),val(2),val(3));       
 end
 ligne=fgetl(fid);  
 iter=iter+1;
 end
end
if OutputFile==1
 fprintf(fidWrite,'0 0. 0. 0.\n');      
elseif OutputFile==2
 fprintf(fidWrite,'*RETURN\n');
 elseif OutputFile==3
 fprintf(fidWrite,'ENDCOORDINATES\n'); 
 fprintf(fidWrite,'\n');
 fprintf(fidWrite,'PANELS TYPE 0\n');    
end
 ligne=fgetl(fid); 
 ligne=fgetl(fid);  
 val=fscanf(fid,'%g',4);
 NElements=val(2);
 ligne=fgetl(fid);
 iter=1;NiterElem=1;
 ibdyline=1;
 while iter<=NElements
     val0=fscanf(fid,'%i',4);
     ligne=fgetl(fid);%
     Neltemp=val0(4);
     for kk=1:Neltemp
         if val0(1)==0              %Points
             val=fscanf(fid,'%i',2);
             ligne=fgetl(fid);
         elseif val0(1)==1          %Lines
             val=fscanf(fid,'%i',3);
             ligne=fgetl(fid);
              if (OutputFile==1 && IncludeBdyLine==1)
                  bdyline(ibdyline,1:2)=[val(2),val(3)];
                  ibdyline=ibdyline+1;
              end
         elseif val0(1)==2          %quadrilateral panels
             val=fscanf(fid,'%i',5);
             ligne=fgetl(fid);
                 if OutputFile==1 
                     fprintf(fidWrite,'%i %i %i %i \n',val(2),val(3),val(4),val(5));
                 else
                     fprintf(fidWrite,'%i %i %i %i \n',val(2)+NPstarts,val(3)+NPstarts,val(4)+NPstarts,val(5)+NPstarts);
                 end
                 NiterElem=NiterElem+1;
         end
         iter=iter+1;
     end
 end
 fclose(fid);
  if (OutputFile==1 && IncludeBdyLine==1)
     for ii=1:ibdyline-1
       fprintf(fidWrite,'%i %i\n', bdyline(ii,1),bdyline(ii,2));  
     end
  end
if OutputFile==1
    fidW=fopen([dir,filesep,MeshfileOutname,'_FSNem_info.dat'],'w');
    fprintf(fidW,'Number of points    = %i \n',NNodes);
    fprintf(fidW,'Number of panels    = %i \n',NiterElem-1);
    if  IncludeBdyLine==1
    fprintf(fidW,'Number of bdy lines = %i \n',ibdyline-1);
    end
    fclose(fidW);
 fprintf(fidWrite,'0 0. 0. 0.\n');      
elseif OutputFile==2
 fprintf(fidWrite,'*RETURN\n');
 
 fidW=fopen([dir,filesep,MeshfileOutname,'_FSHS_info.dat'],'w');
 fprintf(fidW,'Number of points    = %i \n',NNodes);
 fprintf(fidW,'Number of panels    = %i \n',NiterElem-1);
 if  IncludeBdyLine==1
     fprintf(fidW,'Number of bdy lines = %i \n',ibdyline-1);
 end
 fclose(fidW);
 
 elseif OutputFile==3
 fprintf(fidWrite,'ENDPANELS\n');   
 fprintf(fidWrite,'\n'); 
 fprintf(fidWrite,'ENDFILE\n');
 end
fclose(fidWrite);

disp(['Number of nodes     = ',num2str(NNodes)]);
disp(['Number of panels    = ',num2str(NiterElem-1)]);
disp(['Number of bdy lines = ',num2str(ibdyline-1)]);

disp('done!')
if OutputFile==2
disp('Copy the mesh data to inputmesh.dat file for the hydrostar. Combined with the body mesh data') 
end

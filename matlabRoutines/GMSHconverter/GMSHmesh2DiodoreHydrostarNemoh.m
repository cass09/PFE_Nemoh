%Mesh converter from GMSH to NEMOH, DIODORE OR HYDROSTAR
%Created by Ruddy Kurnia, ECN-LHEAA 2021-2022
%Note only quadrilateral mesh allow

clear all;clc;
OutputFile=1;%1 for Nemoh, 2 for Diodore, 3 for Hydrostar
Meshfile='OC3_Hywind_978'; % input file name without ext
MeshfileOutname='OC3_Hywind_978'; % output file name without ext
RemovedLid=0; %1, if Any it will be removed. If OutputFile=2 (Diodore), Lid always be removed
dir='C:\Users\rkurnia\Documents\Test-VirtualBox\gmsh-4.8.4-Linux64\bin\OC3_Hywind';
fid=fopen([dir,filesep,Meshfile,'.msh'],'r');
if OutputFile==1
IdSym=input('\n -input 1 if half-body mesh, 0 otherwise: ');
fidWrite=fopen([dir,filesep,MeshfileOutname,'_Nem.dat'],'w');
fprintf(fidWrite,'%i %i\n',2,IdSym); %%  0 is not symmetric, 1 is symmetric  
elseif OutputFile==2
fidWrite=fopen([dir,filesep,MeshfileOutname,'_Dio.DAT'],'w');
elseif OutputFile==3
fidWrite=fopen([dir,filesep,MeshfileOutname,'_HS.dat'],'w');
fprintf(fidWrite,'PROJECT xxxxxxxxxxxxxx \n');
fprintf(fidWrite,'USER    xxxxxxxxxxxxxx \n');
fprintf(fidWrite,'\n');

fprintf(fidWrite,'NBBODY  1\n');
fprintf(fidWrite,'SYMMETRY_BODY     1 0\n');
fprintf(fidWrite,'REFLENGTH         1.0\n');
fprintf(fidWrite,'RHO               1025.0\n');
fprintf(fidWrite,'GRAVITY           9.81\n');
fprintf(fidWrite,'\n');
fprintf(fidWrite,'COORDINATES TYPE 0\n');
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
 fprintf(fidWrite,'%i %g %g	%g\n',iter,val(1),val(2),val(3));   
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
 while iter<=NElements
     val0=fscanf(fid,'%i',4);
     ligne=fgetl(fid);%
     Neltemp=val0(4);
     for kk=1:Neltemp
         if val0(1)==0
             val=fscanf(fid,'%i',2);
             ligne=fgetl(fid);
         elseif val0(1)==1
             val=fscanf(fid,'%i',3);
             ligne=fgetl(fid);
         elseif val0(1)==2
             val=fscanf(fid,'%i',5);
             ligne=fgetl(fid);
             P=val(2:5);
             Vect1=[XYZ(P(2),1)-XYZ(P(1),1) XYZ(P(2),2)-XYZ(P(1),2) XYZ(P(2),3)-XYZ(P(1),3)];
             Vect2=[XYZ(P(4),1)-XYZ(P(2),1) XYZ(P(4),2)-XYZ(P(2),2) XYZ(P(4),3)-XYZ(P(2),3)];
             n1=cross(Vect1,Vect2);A1=sqrt(n1(1)^2+n1(2)^2+n1(3)^2)/2;
             Vect3=[XYZ(P(4),1)-XYZ(P(3),1) XYZ(P(4),2)-XYZ(P(3),2) XYZ(P(4),3)-XYZ(P(3),3)];
             Vect4=[XYZ(P(2),1)-XYZ(P(3),1) XYZ(P(2),2)-XYZ(P(3),2) XYZ(P(2),3)-XYZ(P(3),3)];
             n2=cross(Vect3,Vect4);A2=sqrt(n2(1)^2+n2(2)^2+n2(3)^2)/2;
             n=n1+n2;norm_n=sqrt(n(1)^2+n(2)^2+n(3)^2);
             n=n/norm_n;
             Area=A1+A2;
             XM=1/3*(XYZ(P(1),1)+XYZ(P(2),1)+XYZ(P(4),1))*A1/Area+1/3*(XYZ(P(2),1)+XYZ(P(3),1)+XYZ(P(4),1))*A2/Area;
             YM=1/3*(XYZ(P(1),2)+XYZ(P(2),2)+XYZ(P(4),2))*A1/Area+1/3*(XYZ(P(2),2)+XYZ(P(3),2)+XYZ(P(4),2))*A2/Area;
             ZM=1/3*(XYZ(P(1),3)+XYZ(P(2),3)+XYZ(P(4),3))*A1/Area+1/3*(XYZ(P(2),3)+XYZ(P(3),3)+XYZ(P(4),3))*A2/Area;
             if ZM<0 % if not lid meshes
                 if OutputFile==1 || OutputFile==3
                     fprintf(fidWrite,'%i %i %i %i \n',val(2),val(3),val(4),val(5));
                 elseif OutputFile==2
                     fprintf(fidWrite,'%i %i %i	%i %i \n',val(1),val(2),val(3),val(4),val(5));
                 end
                 NiterElem=NiterElem+1;
             elseif ZM==0&& (OutputFile==1 && RemovedLid==0)
                     fprintf(fidWrite,'%i %i %i %i \n',val(2),val(3),val(4),val(5));
                      NiterElem=NiterElem+1;
             elseif ZM>0 && OutputFile==2
                     fprintf(fidWrite,'%i %i %i	%i %i \n',val(1),val(2),val(3),val(4),val(5));
                      NiterElem=NiterElem+1;
             end
         end
         iter=iter+1;
     end
 end
 fclose(fid);
if OutputFile==1
    fidW=fopen([dir,filesep,MeshfileOutname,'_Nem_info.dat'],'w');
    fprintf(fidW,'%i %i  Number of points and number of panels\n',NNodes,NiterElem-1);
    fclose(fidW);

 fprintf(fidWrite,'0 0. 0. 0.\n');      
elseif OutputFile==2
 fprintf(fidWrite,'*RETURN\n');
 elseif OutputFile==3
 fprintf(fidWrite,'ENDPANELS\n');   
 fprintf(fidWrite,'\n'); 
 fprintf(fidWrite,'ENDFILE\n');
 fidW=fopen([dir,filesep,MeshfileOutname,'_HS_info.dat'],'w');
 fprintf(fidW,'%i %i  Number of points and number of panels\n',NNodes,NiterElem-1);
 fclose(fidW);
end
fclose(fidWrite);
disp([num2str(NNodes),' nodes and ',num2str(NiterElem-1),' panels.'])
disp('done!')

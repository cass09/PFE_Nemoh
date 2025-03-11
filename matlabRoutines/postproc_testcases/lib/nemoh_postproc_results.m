clear all;close all
axesWT=2; %0=omega [rad/s], 1=frequency [1/s], 2=period [s];
ExctFdat=load('./diodore_results/excitation_forces.txt','-ascii');
ExctFphasedat=load('./diodore_results/excitation_forces_phase.txt','-ascii');

Aj1dat=load('./diodore_results/addedmass_line1.txt','-ascii');
Aj3dat=load('./diodore_results/addedmass_line3.txt','-ascii');
Aj5dat=load('./diodore_results/addedmass_line5.txt','-ascii');
Bj1dat=load('./diodore_results/damping_line1.txt','-ascii');
Bj3dat=load('./diodore_results/damping_line3.txt','-ascii');
Bj5dat=load('./diodore_results/damping_line5.txt','-ascii');

Idlimitaxes=0;
xlimval=[0.054 0.212];

%IddeflM=5;% if there is deflection mode please specify the mode ID

fid=fopen('ID.dat');
line=fgetl(fid);
rep=fscanf(fid,'%s',1);
fclose('all');
fid=fopen([rep,filesep,'Nemoh.cal'],'r');
for i=1:6
    ligne=fgetl(fid);
end
nBodies=fscanf(fid,'%g',1);
 ligne=fgetl(fid);
for i=1:nBodies
    for j=1:3
        ligne=fgetl(fid);
    end
    Ndof=fscanf(fid,'%g',1);
    ligne=fgetl(fid);
    for j=1:Ndof
    ligne=fgetl(fid);    
    end
    Nforce=fscanf(fid,'%g',1);
    ligne=fgetl(fid);
    for j=1:Nforce
    ligne=fgetl(fid);    
    end
    Nadd=fscanf(fid,'%g',1);
    ligne=fgetl(fid);
    if Nadd>0
        Ndeflm=fscanf(fid,'%g',1);
        ligne=fgetl(fid);
        for j=1:Nadd-1
           ligne=fgetl(fid);
        end
    else
        Ndeflm=0;
    end    
end
ligne=fgetl(fid);
nw=fscanf(fid,'%g',1);
fclose(fid);
fid=fopen([rep,filesep,'results',filesep,'ExcitationForce.tec'],'r');
ligne=fgetl(fid);
Nforce=Nforce+Ndeflm;
Ndof=Ndof+Ndeflm;
for c=1:Nforce*nBodies
    ligne=fgetl(fid);
end
 ligne=fgetl(fid);
for k=1:nw
    ligne=fscanf(fid,'%f',1+2*Nforce*nBodies);
    w(k)=ligne(1);
    for j=1:Nforce*nBodies
        Famp(k,j)=ligne(2*j);
        Fphi(k,j)=ligne(2*j+1);
    end
end
status=fclose(fid);
fid=fopen([rep,filesep,'results',filesep,'RadiationCoefficients.tec'],'r');
ligne=fgetl(fid);
for i=1:Ndof*nBodies
    ligne=fgetl(fid);
end
for i=1:nBodies*Ndof
    ligne=fgetl(fid);
    for k=1:nw
        ligne=fscanf(fid,'%f',1+2*Ndof*nBodies);
        for j=1:Ndof*nBodies
            A(i,j,k)=ligne(2*j);
            B(i,j,k)=ligne(2*j+1);
        end
        ligne=fgetl(fid);
    end
end
status=fclose(fid);
% Expression des efforts d excitation de houle sous la forme Famp*exp(i*Fphi)
Fe=Famp.*(cos(Fphi)+1i*sin(Fphi));
% 
if axesWT==0
    axw=w;
    axwdiod=2*pi./Aj1dat(:,1);
    xlab='\omega [rad/s]';
elseif axesWT==1
    axwdiod=1./Aj1dat(:,1);
    axw=w/2/pi;
    xlab='f [1/s]';    
else
    axwdiod=Aj1dat(:,1);
    axw=2*pi./w;
    xlab='T[s]';
end
a11=squeeze(squeeze(A(1,1,:)));a13=squeeze(squeeze(A(1,3,:)));a15=squeeze(squeeze(A(1,5,:)));
a31=squeeze(squeeze(A(3,1,:)));a33=squeeze(squeeze(A(3,3,:)));a35=squeeze(squeeze(A(3,5,:)));
a51=squeeze(squeeze(A(5,1,:)));a53=squeeze(squeeze(A(5,3,:)));a55=squeeze(squeeze(A(5,5,:)));
b11=squeeze(squeeze(B(1,1,:)));b13=squeeze(squeeze(B(1,3,:)));b15=squeeze(squeeze(B(1,5,:)));
b31=squeeze(squeeze(B(3,1,:)));b33=squeeze(squeeze(B(3,3,:)));b35=squeeze(squeeze(B(3,5,:)));
b51=squeeze(squeeze(B(5,1,:)));b53=squeeze(squeeze(B(5,3,:)));b55=squeeze(squeeze(B(5,5,:)));

figure
subplot(3,3,1)
plot(axwdiod,Aj1dat(:,2),'b-',axw,a11,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel('a_{11}' )
plot_properties;
subplot(3,3,2)
plot(axwdiod,Aj1dat(:,4),'b-',axw,a13,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel('a_{13}' )
plot_properties;
subplot(3,3,3)
plot(axwdiod,Aj1dat(:,6),'b-',axw,a15,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel('a_{15}' )
plot_properties;
subplot(3,3,4)
plot(axwdiod,Aj3dat(:,2),'b-',axw,a31,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel('a_{31}' )
plot_properties;
subplot(3,3,5)
plot(axwdiod,Aj3dat(:,4),'b-',axw,a33,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel('a_{33}' )
plot_properties;
subplot(3,3,6)
plot(axwdiod,Aj3dat(:,6),'b-',axw,a35,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel('a_{35}' )
plot_properties;
subplot(3,3,7)
plot(axwdiod,Aj5dat(:,2),'b-',axw,a51,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
xlabel(xlab)
grid ON
ylabel('a_{51}' )
plot_properties;
subplot(3,3,8)
plot(axwdiod,Aj5dat(:,4),'b-',axw,a53,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
xlabel(xlab)
grid ON
ylabel('a_{53}' )
plot_properties;
subplot(3,3,9)
plot(axwdiod,Aj5dat(:,6),'b-',axw,a55,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
xlabel(xlab)
grid ON
ylabel('a_{55}' )
plot_properties;

figure
subplot(3,3,1)
plot(axwdiod,Bj1dat(:,2),'b-',axw,b11,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel('b_{11}' )
plot_properties;
subplot(3,3,2)
plot(axwdiod,Bj1dat(:,4),'b-',axw,b13,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel('b_{13}' )
plot_properties;
subplot(3,3,3)
plot(axwdiod,Bj1dat(:,6),'b-',axw,b15,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel('b_{15}' )
plot_properties;
subplot(3,3,4)
plot(axwdiod,Bj3dat(:,2),'b-',axw,b31,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel('b_{31}' )
plot_properties;
subplot(3,3,5)
plot(axwdiod,Bj3dat(:,4),'b-',axw,b33,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel('b_{33}' )
plot_properties;
subplot(3,3,6)
plot(axwdiod,Bj3dat(:,6),'b-',axw,b35,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel('b_{35}' )
plot_properties;
subplot(3,3,7)
plot(axwdiod,Bj5dat(:,2),'b-',axw,b51,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
xlabel(xlab)
grid ON
ylabel('b_{51}' )
plot_properties;
subplot(3,3,8)
plot(axwdiod,Bj5dat(:,4),'b-',axw,b53,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
xlabel(xlab)
grid ON
ylabel('b_{53}' )
plot_properties;
subplot(3,3,9)
plot(axwdiod,Bj5dat(:,6),'b-',axw,b55,'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
xlabel(xlab)
grid ON
ylabel('b_{55}' )
plot_properties;

figure
subplot(3,2,1)
plot(axwdiod,ExctFdat(:,2),'b-',axw,abs(Fe(:,1)),'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel('|F_1|(N)' )
plot_properties;
subplot(3,2,3)
plot(axwdiod,ExctFdat(:,4),'b-',axw,abs(Fe(:,3)),'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel('|F_3|(N)' )
plot_properties;
subplot(3,2,5)
plot(axwdiod,ExctFdat(:,6),'b-',axw,abs(Fe(:,5)),'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel('|F_5|(N)' )
xlabel(xlab)
plot_properties;

subplot(3,2,2)
plot(axwdiod,ExctFphasedat(:,2)+pi/2,'b-',axw,Fphi(:,1),'r--')
legend('Diodore \angle+pi/2','Nemoh')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel(' \angle F_1(rad)' )
plot_properties;
subplot(3,2,4)
plot(axwdiod,ExctFphasedat(:,4)+pi/2,'b-',axw,Fphi(:,3),'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel(' \angle F_3(rad)' )
plot_properties;
subplot(3,2,6)
plot(axwdiod,ExctFphasedat(:,6)+pi/2,'b-',axw,Fphi(:,5),'r--')
if Idlimitaxes==1
   xlim([xlimval]); 
end
%xlabel('T (s)')
grid ON
ylabel(' \angle F_5(rad)' )
xlabel(xlab)
plot_properties;



% if Ndeflm>0
%    %%Plotting
%    Ndof=Ndof-Ndeflm;
% 
%    a1m=squeeze(squeeze(A(1,Ndof+IddeflM,:)));a3m=squeeze(squeeze(A(3,Ndof+IddeflM,:)));a5m=squeeze(squeeze(A(5,Ndof+IddeflM,:)));amm=squeeze(squeeze(A(Ndof+IddeflM,Ndof+IddeflM,:)));
%    a1m1=squeeze(squeeze(A(1,Ndof+IddeflM+1,:)));a3m1=squeeze(squeeze(A(3,Ndof+IddeflM+1,:)));a5m1=squeeze(squeeze(A(5,Ndof+IddeflM+1,:)));am1m1=squeeze(squeeze(A(Ndof+IddeflM+1,Ndof+IddeflM+1,:)));
%    a1m2=squeeze(squeeze(A(1,Ndof+IddeflM+2,:)));a3m2=squeeze(squeeze(A(3,Ndof+IddeflM+2,:)));a5m2=squeeze(squeeze(A(5,Ndof+IddeflM+2,:)));am2m2=squeeze(squeeze(A(Ndof+IddeflM+2,Ndof+IddeflM+2,:)));
%    a1m3=squeeze(squeeze(A(1,Ndof+IddeflM+3,:)));a3m3=squeeze(squeeze(A(3,Ndof+IddeflM+3,:)));a5m3=squeeze(squeeze(A(5,Ndof+IddeflM+3,:)));am3m3=squeeze(squeeze(A(Ndof+IddeflM+3,Ndof+IddeflM+3,:)));
%   
%    b1m=squeeze(squeeze(B(1,Ndof+IddeflM,:)));b3m=squeeze(squeeze(B(3,Ndof+IddeflM,:)));b5m=squeeze(squeeze(B(5,Ndof+IddeflM,:)));bmm=squeeze(squeeze(B(Ndof+IddeflM,Ndof+IddeflM,:)));
%    b1m1=squeeze(squeeze(B(1,Ndof+IddeflM+1,:)));b3m1=squeeze(squeeze(B(3,Ndof+IddeflM+1,:)));b5m1=squeeze(squeeze(B(5,Ndof+IddeflM+1,:)));bm1m1=squeeze(squeeze(B(Ndof+IddeflM+1,Ndof+IddeflM+1,:)));
%    b1m2=squeeze(squeeze(B(1,Ndof+IddeflM+2,:)));b3m2=squeeze(squeeze(B(3,Ndof+IddeflM+2,:)));b5m2=squeeze(squeeze(B(5,Ndof+IddeflM+2,:)));bm2m2=squeeze(squeeze(B(Ndof+IddeflM+2,Ndof+IddeflM+2,:)));
%    b1m3=squeeze(squeeze(B(1,Ndof+IddeflM+3,:)));b3m3=squeeze(squeeze(B(3,Ndof+IddeflM+3,:)));b5m3=squeeze(squeeze(B(5,Ndof+IddeflM+3,:)));bm3m3=squeeze(squeeze(B(Ndof+IddeflM+3,Ndof+IddeflM+3,:)));
% 
% 
% figure
% subplot(4,4,1)
% plot(axw,a1m,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['a_{1,',num2str(Ndof+IddeflM),'}'])
% plot_properties;
% subplot(4,4,2)
% plot(axw,a3m,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['a_{3,',num2str(Ndof+IddeflM),'}'])
% plot_properties;
% subplot(4,4,3)
% plot(axw,a5m,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['a_{5,',num2str(Ndof+IddeflM),'}'])
% plot_properties;
% subplot(4,4,4)
% plot(axw,amm,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['a_{',num2str(Ndof+IddeflM),',',num2str(Ndof+IddeflM),'}'])
% plot_properties;
% subplot(4,4,5)
% plot(axw,a1m1,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['a_{1,',num2str(Ndof+IddeflM+1),'}'])
% plot_properties;
% subplot(4,4,6)
% plot(axw,a3m1,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['a_{3,',num2str(Ndof+IddeflM+1),'}'])
% plot_properties;
% subplot(4,4,7)
% plot(axw,a5m1,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['a_{5,',num2str(Ndof+IddeflM+1),'}'])
% plot_properties;
% subplot(4,4,8)
% plot(axw,am1m1,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['a_{',num2str(Ndof+IddeflM+1),',',num2str(Ndof+IddeflM+1),'}'])
% plot_properties;
% subplot(4,4,9)
% plot(axw,a1m2,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['a_{1,',num2str(Ndof+IddeflM+2),'}'])
% plot_properties;
% subplot(4,4,10)
% plot(axw,a3m2,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['a_{3,',num2str(Ndof+IddeflM+2),'}'])
% plot_properties;
% subplot(4,4,11)
% plot(axw,a5m2,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['a_{5,',num2str(Ndof+IddeflM+2),'}'])
% plot_properties;
% subplot(4,4,12)
% plot(axw,am2m2,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['a_{',num2str(Ndof+IddeflM+2),',',num2str(Ndof+IddeflM+2),'}'])
% plot_properties;
% subplot(4,4,13)
% plot(axw,a1m3,'r--')
% %xlabel('T (s)')
% xlabel(xlab)
% grid ON
% ylabel(['a_{1,',num2str(Ndof+IddeflM+3),'}'])
% plot_properties;
% subplot(4,4,14)
% plot(axw,a3m3,'r--')
% %xlabel('T (s)')
% xlabel(xlab)
% grid ON
% ylabel(['a_{3,',num2str(Ndof+IddeflM+3),'}'])
% plot_properties;
% subplot(4,4,15)
% plot(axw,a5m3,'r--')
% %xlabel('T (s)')
% xlabel(xlab)
% grid ON
% ylabel(['a_{5,',num2str(Ndof+IddeflM+3),'}'])
% plot_properties;
% subplot(4,4,16)
% plot(axw,am3m3,'r--')
% %xlabel('T (s)')
% xlabel(xlab)
% grid ON
% ylabel(['a_{',num2str(Ndof+IddeflM+3),',',num2str(Ndof+IddeflM+3),'}'])
% plot_properties;
% 
% figure
% subplot(4,4,1)
% plot(axw,b1m,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['b_{1,',num2str(Ndof+IddeflM),'}'])
% plot_properties;
% subplot(4,4,2)
% plot(axw,b3m,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['b_{3,',num2str(Ndof+IddeflM),'}'])
% plot_properties;
% subplot(4,4,3)
% plot(axw,b5m,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['b_{5,',num2str(Ndof+IddeflM),'}'])
% plot_properties;
% subplot(4,4,4)
% plot(axw,bmm,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['b_{',num2str(Ndof+IddeflM),',',num2str(Ndof+IddeflM),'}'])
% plot_properties;
% subplot(4,4,5)
% plot(axw,b1m1,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['b_{1,',num2str(Ndof+IddeflM+1),'}'])
% plot_properties;
% subplot(4,4,6)
% plot(axw,b3m1,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['b_{3,',num2str(Ndof+IddeflM+1),'}'])
% plot_properties;
% subplot(4,4,7)
% plot(axw,b5m1,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['b_{5,',num2str(Ndof+IddeflM+1),'}'])
% plot_properties;
% subplot(4,4,8)
% plot(axw,bm1m1,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['b_{',num2str(Ndof+IddeflM+1),',',num2str(Ndof+IddeflM+1),'}'])
% plot_properties;
% subplot(4,4,9)
% plot(axw,b1m2,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['b_{1,',num2str(Ndof+IddeflM+2),'}'])
% plot_properties;
% subplot(4,4,10)
% plot(axw,b3m2,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['b_{3,',num2str(Ndof+IddeflM+2),'}'])
% plot_properties;
% subplot(4,4,11)
% plot(axw,b5m2,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['b_{5,',num2str(Ndof+IddeflM+2),'}'])
% plot_properties;
% subplot(4,4,12)
% plot(axw,bm2m2,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['b_{',num2str(Ndof+IddeflM+2),',',num2str(Ndof+IddeflM+2),'}'])
% plot_properties;
% subplot(4,4,13)
% plot(axw,b1m3,'r--')
% %xlabel('T (s)')
% xlabel(xlab)
% grid ON
% ylabel(['b_{1,',num2str(Ndof+IddeflM+3),'}'])
% xlabel(xlab)
% plot_properties;
% subplot(4,4,14)
% plot(axw,b3m3,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['b_{3,',num2str(Ndof+IddeflM+3),'}'])
% xlabel(xlab)
% plot_properties;
% subplot(4,4,15)
% plot(axw,b5m3,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['b_{5,',num2str(Ndof+IddeflM+3),'}'])
% xlabel(xlab)
% plot_properties;
% subplot(4,4,16)
% plot(axw,bm3m3,'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['b_{',num2str(Ndof+IddeflM+3),',',num2str(Ndof+IddeflM+3),'}'])
% xlabel(xlab)
% plot_properties;
% 
% figure
% subplot(4,1,1)
% plot(axw,abs(Fe(:,Ndof+IddeflM)),'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['F_{',num2str(Ndof+IddeflM),'} (N)'])
% plot_properties;
% subplot(4,1,2)
% plot(axw,abs(Fe(:,Ndof+IddeflM+1)),'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['F_{',num2str(Ndof+IddeflM+1),'} (N)'])
% plot_properties;
% subplot(4,1,3)
% plot(axw,abs(Fe(:,Ndof+IddeflM+2)),'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['F_{',num2str(Ndof+IddeflM+2),'} (N)'])
% plot_properties;
% subplot(4,1,4)
% plot(axw,abs(Fe(:,Ndof+IddeflM+3)),'r--')
% %xlabel('T (s)')
% grid ON
% ylabel(['F_{',num2str(Ndof+IddeflM+3),'} (N)'])
% xlabel(xlab)
% plot_properties;
% end
%%Plotting
% figure
% subplot(3,2,1)
% a(1,:)=A(1,1,:);
% b(1,:)=B(1,1,:);
% plot(axw,a,'b-+',axw,b,'r-+')
% %xlabel('T (s)')
% grid ON
% ylabel('Coeffs' )
% legend('Added mass (A)','Damping (B)')
% title('Surge')
% plot_properties;
% 
% subplot(3,2,2)
% plot(axw,abs(Fe(:,1)),'g-+')
% %xlabel('T (s)')
% grid ON
% ylabel('F (N)' )
% title('Surge')
% plot_properties;
% 
% subplot(3,2,3)
% a(1,:)=A(3,3,:);
% b(1,:)=B(3,3,:);
% plot(axw,a,'b-+',axw,b,'r-+')
% %xlabel('T (s)')
% grid ON
% ylabel('Coeffs' )
% %legend('Added mass (A)','Damping (B)')
% title('Heave')
% plot_properties;
% 
% subplot(3,2,4)
% plot(axw,abs(Fe(:,3)),'g-+')
% %xlabel('T (s)')
% grid ON
% ylabel('F(N)' )
% title('Heave')
% plot_properties;
% 
% subplot(3,2,5)
% a(1,:)=A(5,5,:);
% b(1,:)=B(5,5,:);
% plot(axw,a,'b-+',axw,b,'r-+')
% xlabel(xlab)
% ylabel('Coeffs' )
% grid ON
% %legend('Added mass (A)','Damping (B)')
% title('Pitch')
% plot_properties;
% 
% subplot(3,2,6)
% plot(axw,abs(Fe(:,5)),'g-+')
% xlabel(xlab)
% grid ON
% ylabel('F(N)' )
% title('Pitch')
% %suptitle('Added Masses and Excitation Forces')
% plot_properties;
% 
% if Ndeflm>0
%    %%Plotting
%    Ndof=Ndof-Ndeflm;
% figure
% subplot(3,2,1)
% a(1,:)=A(Ndof+IddeflM,Ndof+IddeflM,:);
% b(1,:)=B(Ndof+IddeflM,Ndof+IddeflM,:);
% plot(axw,a,'b-+',axw,b,'r-+')
% %xlabel('T (s)')
% grid ON
% ylabel('Coeffs' )
% legend('Added mass (A)','Damping (B)')
% title(['Defl mode ', num2str(IddeflM)])
% plot_properties;
% 
% subplot(3,2,2)
% plot(axw,abs(Fe(:,Ndof+IddeflM)),'g-+')
% %xlabel('T (s)')
% grid ON
% ylabel('F (N)' )
% title(['Defl mode ', num2str(IddeflM)])
% plot_properties;
% 
% subplot(3,2,3)
% a(1,:)=A(Ndof+IddeflM+1,Ndof+IddeflM+1,:);
% b(1,:)=B(Ndof+IddeflM+1,Ndof+IddeflM+1,:);
% plot(axw,a,'b-+',axw,b,'r-+')
% %xlabel('T (s)')
% grid ON
% ylabel('Coeffs' )
% %legend('Added mass (A)','Damping (B)')
% title(['Defl mode ', num2str(IddeflM+1)])
% plot_properties;
% 
% subplot(3,2,4)
% plot(axw,abs(Fe(:,Ndof+IddeflM+1)),'g-+')
% %xlabel('T (s)')
% grid ON
% ylabel('F(N)' )
% title(['Defl mode ', num2str(IddeflM+1)])
% plot_properties;
% 
% subplot(3,2,5)
% a(1,:)=A(Ndof+IddeflM+2,Ndof+IddeflM+2,:);
% b(1,:)=B(Ndof+IddeflM+2,Ndof+IddeflM+2,:);
% plot(axw,a,'b-+',axw,b,'r-+')
% xlabel(xlab)
% ylabel('Coeffs' )
% grid ON
% %legend('Added mass (A)','Damping (B)')
% title(['Defl mode ', num2str(IddeflM+2)])
% plot_properties;
% 
% subplot(3,2,6)
% plot(axw,abs(Fe(:,Ndof+IddeflM+2)),'g-+')
% xlabel(xlab)
% grid ON
% ylabel('F(N)' )
% title(['Defl mode ', num2str(IddeflM+2)])
% %suptitle('Added Masses and Excitation Forces') 
% plot_properties;
% end
%%
% %=========== MOTIONS RAOS FROM CUMMINS EQUATION IN FREQUENCY DOMAIN =====%
% %--- Calculation--------------%
% %close all
% %clear all
% % load Nemohresult_test
% % load Mesh_outputs
% 
% for k=1:length(w)
%     Fe_omega(1,:)=Fe(k,:); A_omega(:,:)= A(:,:,k)  ;B_omega(:,:)= B(:,:,k)  ;
%     RAO(k,:)=(Fe_omega)/(-(M+A_omega)*(w(k))^2-1i*w(k)*(B_omega)+KHyd); % Coupling between DoF
%     for p=1:6
%         RAO1(k,p)=(Fe(k,p))/(-(M(p,p)+A(p,p,k))*(w(k))^2-1i*w(k)*(B(p,p,k))+KHyd(p,p)); % No coupling between the DoF
%     end
% end
% 
% %---------Plot the results--------------%
% freq=w/(2*pi);
% title_RAO ={'Surge' 'Sway' 'Heave' 'Roll' 'Pitch' 'Yaw'};
% figure
% for jj=1:6
% subplot(2,3,jj) % Surge
% plot(freq,abs(RAO(:,jj)),'b-+','LineWidth',2.0)
% title([title_RAO(jj)]);
% xlabel('Frequency(Hz)')
% grid ON
% ylabel('RAO' )
% hold off
% end
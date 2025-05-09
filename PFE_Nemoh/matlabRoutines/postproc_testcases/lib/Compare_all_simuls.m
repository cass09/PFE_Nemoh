clear all;close all
[pathstr,~,~] = fileparts(mfilename('fullpath'));
cd (pathstr);
addpath(genpath(pathstr))

axesWT=0; %0=omega [rad/s], 1=frequency [1/s], 2=period [s], 3=w^2/g
IDDOF=5; % ID dof 1,3,5 for surge, heave, pitch
IDDIOD=0; % if to be compare with diodore 1
IDHYDROSTAR=1;
IDRAO=1;
Idlimitaxes=1;
xlimval=[0 4];
%xlimval=[0.22 0.35];
IdNormalized=0; %% normalized only excitation force by rho g, phase by -pi to be compare with wamit
rho=1025;
g=9.81;

data_file={'SOFTWIND_gmshmesh_half_1872-QTF-COG_redo_eps001-dir10'};
hydrostarDir='hydrostar_results';

lintype={'--r','-.g',':c',':m'};
HySmb='b-';
Ndata_file=length(data_file);
if IDRAO==1
figure('unit','normalized','Position',[0.3 0.1 0.2 0.2]);axRAO=axes;
figure('unit','normalized','Position',[0.3 0.1 0.2 0.2]);axRAOP=axes;
end
figure('unit','normalized','Position',[0.3 0.1 0.2 0.5])
for II=1:Ndata_file

    if IDHYDROSTAR==1
    NHeader=22;
    F1datH=fun_readHydrostarOutput(['./',hydrostarDir,'/fxf1st.rao'],NHeader);
    F3datH=fun_readHydrostarOutput(['./',hydrostarDir,'/fzf1st.rao'],NHeader);
    F5datH=fun_readHydrostarOutput(['./',hydrostarDir,'/myf1st.rao'],NHeader);

    A11datH=fun_readHydrostarOutput(['./',hydrostarDir,'/AddedMass_11.dat'],NHeader);
    A33datH=fun_readHydrostarOutput(['./',hydrostarDir,'/AddedMass_33.dat'],NHeader);
    A55datH=fun_readHydrostarOutput(['./',hydrostarDir,'/AddedMass_55.dat'],NHeader);

    B11datH=fun_readHydrostarOutput(['./',hydrostarDir,'/Damping_11.dat'],NHeader);
    B33datH=fun_readHydrostarOutput(['./',hydrostarDir,'/Damping_33.dat'],NHeader);
    B55datH=fun_readHydrostarOutput(['./',hydrostarDir,'/Damping_55.dat'],NHeader);

%     ExctFdat=load('./diodore_results/excitation_forces.txt','-ascii');
%     ExctFphasedat=load('./diodore_results/excitation_forces_phase.txt','-ascii');

%     Aj1dat=load('./diodore_results/addedmass_line1.txt','-ascii');
%     Aj3dat=load('./diodore_results/addedmass_line3.txt','-ascii');
%     Aj5dat=load('./diodore_results/addedmass_line5.txt','-ascii');
%     Bj1dat=load('./diodore_results/damping_line1.txt','-ascii');
%     Bj3dat=load('./diodore_results/damping_line3.txt','-ascii');
%     Bj5dat=load('./diodore_results/damping_line5.txt','-ascii');
%
%     F1datD=ExctFdat(:,2);F3datD=ExctFdat(:,4);F5datD=ExctFdat(:,6);
%     A11datD=Aj1dat(:,2);A33datD=Aj3dat(:,4);A55datD=Aj5dat(:,6);
%     B11datD=Bj1dat(:,2);B33datD=Bj3dat(:,4);B55datD=Bj5dat(:,6);
    end

    dirfile=(data_file{II});
    cd(dirfile);
    rep='./';

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
    ligne=fscanf(fid,'%d',2);
    nw=ligne(2);
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
        ligne=fgetl(fid);
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
    if IDHYDROSTAR==1
    F1datH(:,3)=deg2rad(F1datH(:,3));
    F3datH(:,3)=deg2rad(F3datH(:,3));
    F5datH(:,3)=deg2rad(F5datH(:,3));
    end
    if IdNormalized==1
        Famp=Famp./rho./g;
        Fphi=wrapToPi(Fphi)./pi;
        if IDHYDROSTAR==1
        F1datH(:,2)=F1datH(:,2)./rho./g;
        F1datH(:,3)=wrapToPi(F1datH(:,3))./pi;
        F3datH(:,2)=F3datH(:,2)./rho./g;
        F3datH(:,3)=wrapToPi(F3datH(:,3))./pi;
        F5datH(:,2)=F5datH(:,2)./rho./g;
        F5datH(:,3)=wrapToPi(F5datH(:,3))./pi;
        end
    end

    if axesWT==0
        axw=w;
        if IDHYDROSTAR==1
        axwH=A11datH(:,1);
        end
        if IDDIOD==1
            axwD=2*pi./Aj1dat(:,1);
        end
        xlab='\omega [rad/s]';
    elseif axesWT==1
        if IDHYDROSTAR==1
        axwH=A11datH(:,1)/2/pi;
        end
        if IDDIOD==1
        axwD=1./Aj1dat(:,1);
        end
        axw=w/2/pi;
        xlab='f [Hz]';
    elseif axesWT==2
        if IDHYDROSTAR==1
        axwH=2*pi./A11datH(:,1);
        end
        if IDDIOD==1
        axwD=Aj1dat(:,1);
        end
        axw=2*pi./w;
        xlab='T[s]';
    elseif axesWT==3
        axw=w.^2/g;
        if IDHYDROSTAR==1
        axwH=A11datH(:,1).^2/g;
        end
        if IDDIOD==1
        axwD=(2*pi./Aj1dat(:,1)).^2/g;
        end
        xlab='\omega^2/g';
    end
    a11=squeeze(squeeze(A(1,1,:)));a13=squeeze(squeeze(A(1,3,:)));a15=squeeze(squeeze(A(1,5,:)));
    a31=squeeze(squeeze(A(3,1,:)));a33=squeeze(squeeze(A(3,3,:)));a35=squeeze(squeeze(A(3,5,:)));
    a51=squeeze(squeeze(A(5,1,:)));a53=squeeze(squeeze(A(5,3,:)));a55=squeeze(squeeze(A(5,5,:)));
    b11=squeeze(squeeze(B(1,1,:)));b13=squeeze(squeeze(B(1,3,:)));b15=squeeze(squeeze(B(1,5,:)));
    b31=squeeze(squeeze(B(3,1,:)));b33=squeeze(squeeze(B(3,3,:)));b35=squeeze(squeeze(B(3,5,:)));
    b51=squeeze(squeeze(B(5,1,:)));b53=squeeze(squeeze(B(5,3,:)));b55=squeeze(squeeze(B(5,5,:)));

    if IDDIOD==1

    else
        subplot(4,1,1)
        if II==1 && IDHYDROSTAR==1
            if IDDOF==1
                plot(axwH,A11datH(:,2),HySmb,axw,a11,lintype{II})
                ylabel('a_{11}' )
            elseif IDDOF==3
                plot(axwH,A33datH(:,2),HySmb,axw,a33,lintype{II})
                ylabel('a_{33}' )
            elseif IDDOF==5
                plot(axwH,A55datH(:,2),HySmb,axw,a55,lintype{II})
                ylabel('a_{55}' )
            end
        else
            if IDDOF==1
                plot(axw,a11,lintype{II})
                ylabel('a_{11}' )
            elseif IDDOF==3
                plot(axw,a33,lintype{II})
                ylabel('a_{33}' )
            elseif IDDOF==5
                plot(axw,a55,lintype{II})
                ylabel('a_{55}' )
            end
        end
        hold on;
        if Idlimitaxes==1
            xlim([xlimval]);
        end
        %xlabel(xlab)
        grid ON

        plot_properties;
        subplot(4,1,2)
        if II==1 && IDHYDROSTAR==1
            if IDDOF==1
                plot(axwH,B11datH(:,2),HySmb,axw,b11,lintype{II})
                ylabel('b_{11}' )
            elseif IDDOF==3
                plot(axwH,B33datH(:,2),HySmb,axw,b33,lintype{II})
                ylabel('b_{33}' )
            elseif IDDOF==5
                plot(axwH,B55datH(:,2),HySmb,axw,b55,lintype{II})
                ylabel('b_{55}' )
            end
        else
             if IDDOF==1
                plot(axw,b11,lintype{II})
                ylabel('b_{11}' )
            elseif IDDOF==3
                plot(axw,b33,lintype{II})
                ylabel('b_{33}' )
            elseif IDDOF==5
                plot(axw,b55,lintype{II})
                ylabel('b_{55}' )
             end
        end
        hold on;
        if Idlimitaxes==1
            xlim([xlimval]);
        end
        %xlabel('T (s)')
        grid ON
        %xlabel(xlab)
        plot_properties;
        subplot(4,1,3)
        if II==1 && IDHYDROSTAR==1
        if IDDOF==1
            plot(axwH,F1datH(:,2),HySmb,axw,Famp(:,1),lintype{II})
            ylabel('|F_{1}|' )
        elseif IDDOF==3
            plot(axwH,F3datH(:,2),HySmb,axw,Famp(:,3),lintype{II})
            ylabel('|F_{3}|' )
        elseif IDDOF==5
            plot(axwH,F5datH(:,2),HySmb,axw,Famp(:,5),lintype{II})
            ylabel('|F_{5}|' )
        end
        else
        if IDDOF==1
            plot(axw,Famp(:,1),lintype{II})
            ylabel('|F_{1}|' )
        elseif IDDOF==3
            plot(axw,Famp(:,3),lintype{II})
            ylabel('|F_{3}|' )
        elseif IDDOF==5
            plot(axw,Famp(:,5),lintype{II})
            ylabel('|F_{5}|' )
        end
        end
        hold on;
        xlabel(xlab)
        if Idlimitaxes==1
            xlim([xlimval]);
        end
        %xlabel('T (s)')

        plot_properties;
        grid ON

        subplot(4,1,4)
        if II==1 && IDHYDROSTAR==1
        if IDDOF==1
            plot(axwH,F1datH(:,3),HySmb,axw,Fphi(:,1),lintype{II})
            ylabel('Phase F_{1}' )
        elseif IDDOF==3
            plot(axwH,F3datH(:,3),HySmb,axw,Fphi(:,3),lintype{II})
            ylabel('Phase F_{3}' )
        elseif IDDOF==5
            plot(axwH,F5datH(:,3),HySmb,axw,Fphi(:,5),lintype{II})
            ylabel('Phase F_{5}' )
        end
        else
        if IDDOF==1
            plot(axw,Fphi(:,1),lintype{II})
            ylabel('Phase F_{1}' )
        elseif IDDOF==3
            plot(axw,Fphi(:,3),lintype{II})
            ylabel('Phase F_{3}' )
        elseif IDDOF==5
            plot(axw,Fphi(:,5),lintype{II})
            ylabel('Phase F_{5}' )
        end
        end
        hold on;
        xlabel(xlab)
        if Idlimitaxes==1
            xlim([xlimval]);
        end
        %xlabel('T (s)')

        plot_properties;
        grid ON
    end
   if IDRAO==1
        RAOAmpl=zeros(nw,6);
        RAOPhase=zeros(nw,6);

        fid=fopen([rep,filesep,'Motion',filesep,'RAO.dat'],'r');

        ligne=fgetl(fid);
        %for ibeta=1:nbeta
            ligne=fgetl(fid);
            for iw=1:nw
                ligne=fscanf(fid,'%g',13);
                RAOAmpl(iw,:)=ligne(2:7);
                RAOPhase(iw,:)=ligne(8:13);
                ligne=fgetl(fid);
            end
        %end
        status=fclose(fid);



        if IDHYDROSTAR==1
        RAOH_surge=fun_readHydrostarOutput(['./../',hydrostarDir,'/surge.rao'],NHeader);
        RAOH_heave=fun_readHydrostarOutput(['./../',hydrostarDir,'/heave.rao'],NHeader);
        RAOH_pitch=fun_readHydrostarOutput(['./../',hydrostarDir,'/pitch.rao'],NHeader);
        end



        if II==1 && IDHYDROSTAR==1
            if IDDOF==1
                plot(axRAO,axwH,RAOH_surge(:,2),HySmb,axw,RAOAmpl(:,1),lintype{II})
                ylabel(axRAO,'|\xi_1/A|' )
                plot(axRAOP,axwH,mod(RAOH_surge(:,3),360),HySmb,axw,RAOPhase(:,1),lintype{II})
                ylabel(axRAOP,'Phase' )
            elseif IDDOF==3
                plot(axRAO,axwH,RAOH_heave(:,2),HySmb,axw,RAOAmpl(:,3),lintype{II})
                ylabel(axRAO,'|\xi_3/A|' )
                plot(axRAOP,axwH,mod(RAOH_heave(:,3),360),HySmb,axw,RAOPhase(:,3),lintype{II})
                ylabel(axRAOP,'Phase' )
            elseif IDDOF==5
                plot(axRAO,axwH,RAOH_pitch(:,2),HySmb,axw,RAOAmpl(:,5),lintype{II})
                ylabel(axRAO,'|\xi_5/A|' )
                plot(axRAOP,axwH,mod(RAOH_pitch(:,3),360),HySmb,axw,RAOPhase(:,5),lintype{II})
                ylabel(axRAOP,'Phase' )
            end
        else
           if IDDOF==1
                plot(axRAO,axw,RAOAmpl(:,1),lintype{II})
                ylabel(axRAO,'|\xi_1/A|' )
                plot(axRAOP,axw,RAOPhase(:,1),lintype{II})
                ylabel(axRAOP,'Phase' )
            elseif IDDOF==3
                plot(axRAO,axw,RAOAmpl(:,3),lintype{II})
                ylabel(axRAO,'|\xi_3/A|' )
                plot(axRAOP,axw,RAOPhase(:,3),lintype{II})
                ylabel(axRAOP,'Phase' )
            elseif IDDOF==5
                plot(axRAO,axw,RAOAmpl(:,5),lintype{II})
                ylabel(axRAO,'|\xi_5/A|' )
                plot(axRAOP,axw,RAOPhase(:,5),lintype{II})
                ylabel(axRAOP,'Phase' )
            end
        end
        % ylim(axRAO,[0,1.5])
         %xlim(axRAO,[0,2])
           xlabel(axRAO,xlab)
             plot_properties;
              hold(axRAO,'on');

             % xlim(axRAOP,[0,2])
           xlabel(axRAOP,xlab)
             plot_properties;
              hold(axRAOP,'on');


    end

    cd ..
end

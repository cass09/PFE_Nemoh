%code for comparing the hydrodynamic coefficients between NEMOH & Hydrostar
%data parameters are set in data_param;
%make sure the testcases data available after run.
%R. Kurnia, 2022 ECN

clear all;
close all
[pathstr,~,~] = fileparts(mfilename('fullpath'));
cd (pathstr);
addpath(genpath(pathstr))

TestCaseNo=1;
IDRAO=0;               % 0 or 1; 0 not computed
if TestCaseNo==1       %3_NonSymmetrical
    dataparamNEMOH1_5;
    IDRAO=0;
elseif TestCaseNo==2   %8a_Cylinder_irregfreq
    dataparamNEMOH1_6;
    IDRAO=0;
elseif  TestCaseNo==3  %8b_QTF_Cylinder
    dataparamNEMOH1_1;
    IDRAO=1;
elseif TestCaseNo==4   %9_QTF_OC4_Semisubmersible
    dataparamNEMOH1_2;
    IDRAO=1;
elseif TestCaseNo==5   %10a_QTF_SOFTWIND
    dataparamNEMOH1_3;
    IDRAO=1;
elseif TestCaseNo==6   %10b_QTF_SOFTWIND_FS
    dataparamNEMOH1_4
    IDRAO=1;
elseif TestCaseNo==7   %11_QTF_OC3_Hywind
    dataparamNEMOH1_7
    IDRAO=1;
end

ID_MODES_PLOT=3;%1 all modes; 2= modes 1,3,5; 3= modes 1,3,5, 2,4,6
axesWT=0; %0=omega [rad/s], 1=frequency [1/s], 2=period [s], 3=w^2/g
IDHYDROSTAR=1;
Idlimitaxes=0;
xlimval=[0 4];%[0 2*pi];
%xlimval=[0.22 0.35];

lintype={'--r','-.g',':c',':m'};
HySmb='b-';
Ndata_file=length(data_file);

if IDHYDROSTAR==1
    for IDDOF=1:6
        if IDDOF==1,    FHdof='fx';
        elseif IDDOF==2,FHdof='fy';
        elseif IDDOF==3,FHdof='fz';
        elseif IDDOF==4,FHdof='mx';
        elseif IDDOF==5,FHdof='my';
        elseif IDDOF==6,FHdof='mz';
        end

        NHeader=22;
        for Idbeta=1:NbetaHdat
            datFH=fun_readMdirHydrostarOutput(...
                ['./',hydrostarDir,'/',FHdof,'f1st.rao'],NHeader,Idbeta,NbetaHdat);

            wwH=datFH(:,1);
            FdatH(:,IDDOF,Idbeta)=datFH(:,2);
            FdatH(:,IDDOF+6,Idbeta)=deg2rad(datFH(:,3));
        end
    end

    AMatH=zeros(size(wwH,1),6,6);
    BMatH=zeros(size(wwH,1),6,6);
    for IDDOF=1:6
        for j=1:IDDOF
            datHA=fun_readMdirHydrostarOutput(['./',hydrostarDir,'/AddedMass_'...
                ,num2str(IDDOF),num2str(j),'.dat'],NHeader,1,NbetaHdat);
            wwH=datHA(:,1);
            AMatH(:,IDDOF,j)=datHA(:,2);
            datHB=fun_readMdirHydrostarOutput(['./',hydrostarDir,'/Damping_'...
                ,num2str(IDDOF),num2str(j),'.dat'],NHeader,1,NbetaHdat);
            BMatH(:,IDDOF,j)=datHB(:,2);
        end
    end


    if IDRAO==1
        for IDDOF=1:6
            if IDDOF==1,    RHdof='surge';
            elseif IDDOF==2,RHdof='sway';
            elseif IDDOF==3,RHdof='heave';
            elseif IDDOF==4,RHdof='roll';
            elseif IDDOF==5,RHdof='pitch';
            elseif IDDOF==6,RHdof='yaw';
            end
            for Idbeta=1:NbetaHdat
                datRAO=fun_readMdirHydrostarOutput(['./',hydrostarDir,'/',RHdof,...
                    '.rao'],NHeader,Idbeta,NbetaHdat);
                wwH=datRAO(:,1);
                RAOH(:,IDDOF,Idbeta)=datRAO(:,2);
                RAOH(:,IDDOF+6,Idbeta)=deg2rad(datRAO(:,3));
            end
        end
    end
end

figure('unit','normalized','Position',[0 0 0.8 0.8]);
for II=1:Ndata_file
    dirfile=[NemohDir,filesep,data_file{II}];
    [wwN,~,AmatNem,~,~]=fun_readHydcoef_NEMOHOutput(dirfile,0,1,0);
    [axw,axwH,xlab]=fun_axes_Adjusted(wwN,wwH,axesWT,IDHYDROSTAR);
    fun_plot_addedMass_dampCoef(axwH,axw,AMatH,AmatNem,HySmb,II,...
        lintype,Idlimitaxes,xlimval,xlab,'a',IDHYDROSTAR,ID_MODES_PLOT)
end

figure('unit','normalized','Position',[0 0 0.8 0.8]);
for II=1:Ndata_file
    dirfile=[NemohDir,filesep,data_file{II}];
    [wwN,~,~,BmatNem,~]=fun_readHydcoef_NEMOHOutput(dirfile,0,1,0);
    [axw,axwH,xlab]=fun_axes_Adjusted(wwN,wwH,axesWT,IDHYDROSTAR);
    fun_plot_addedMass_dampCoef(axwH,axw,BMatH,BmatNem,HySmb,II,...
        lintype,Idlimitaxes,xlimval,xlab,'b',IDHYDROSTAR,ID_MODES_PLOT)
end

figure('unit','normalized','Position',[0 0 0.8 0.8]);
for II=1:Ndata_file
    dirfile=[NemohDir,filesep,data_file{II}];
    [wwN,FdatNem,~,~,~]=fun_readHydcoef_NEMOHOutput(dirfile,1,0,0);
    [axw,axwH,xlab]=fun_axes_Adjusted(wwN,wwH,axesWT,IDHYDROSTAR);
    fun_plot_excitationForcesRAO(axwH,axw,FdatH,FdatNem,HySmb,II,...
        lintype,Idlimitaxes,xlimval,xlab,IDHYDROSTAR,IdbetaPlotH,IdbetaPlot,'|F_{ex}','|',ID_MODES_PLOT)
end
if IDRAO==1
    figure('unit','normalized','Position',[0 0 0.8 0.8]);
    for II=1:Ndata_file
        dirfile=[NemohDir,filesep,data_file{II}];
        [wwN,~,~,~,RAONem]=fun_readHydcoef_NEMOHOutput(dirfile,0,0,1);
        [axw,axwH,xlab]=fun_axes_Adjusted(wwN,wwH,axesWT,IDHYDROSTAR);
        fun_plot_excitationForcesRAO(axwH,axw,RAOH,RAONem,HySmb,II,...
            lintype,Idlimitaxes,xlimval,xlab,IDHYDROSTAR,IdbetaPlotH,IdbetaPlot,'|\xi','/A|',ID_MODES_PLOT)
    end
end

function [axw,axwH,xlab]=fun_axes_Adjusted(wwN,wwH,axesWT,IDHYDROSTAR)
if axesWT==0
    axw=wwN;
    if IDHYDROSTAR==1
        axwH=wwH;
    end
    xlab='\omega [rad/s]';
elseif axesWT==1
    if IDHYDROSTAR==1
        axwH=wwH/2/pi;
    end
    axw=wwN/2/pi;
    xlab='f [Hz]';
elseif axesWT==2
    if IDHYDROSTAR==1
        axwH=2*pi./wwH;
    end
    axw=2*pi./wwN;
    xlab='T[s]';
elseif axesWT==3
    axw=wwN.^2/g;
    if IDHYDROSTAR==1
        axwH=wwH.^2/g;
    end
    xlab='\omega^2/g';
end
end

function fun_plot_addedMass_dampCoef(axwH,axw,AMatH,AmatNem,HySmb,II,...
    lintype,Idlimitaxes,xlimval,xlab,ylab,IDHYDROSTAR,ID_MODES_PLOT)
if ID_MODES_PLOT==1
    subplot(6,6,1)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,1,1))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,1,1))),lintype{II})
      else
        plot(axw,squeeze(squeeze(AmatNem(:,1,1))),lintype{II})
    end
        plot_properties

    ylabel([ylab,'_{11}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(6,6,7)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,2,1))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,2,1))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,2,1))),lintype{II})
    end
         plot_properties

    ylabel([ylab,'_{21}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(6,6,8)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,2,2))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,2,2))),lintype{II})
     else
        plot(axw,squeeze(squeeze(AmatNem(:,2,2))),lintype{II})
    end
        plot_properties

    ylabel([ylab,'_{22}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(6,6,13)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,3,1))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,3,1))),lintype{II})
     else
        plot(axw,squeeze(squeeze(AmatNem(:,3,1))),lintype{II})
    end
    plot_properties
    ylabel([ylab,'_{31}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(6,6,14)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,3,2))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,3,2))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,3,2))),lintype{II})
    end
     plot_properties

    ylabel([ylab,'_{32}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(6,6,15)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,3,3))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,3,3))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,3,3))),lintype{II})
    end
    plot_properties

    ylabel([ylab,'_{33}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(6,6,19)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,4,1))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,4,1))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,4,1))),lintype{II})
    end
      plot_properties

    ylabel([ylab,'_{41}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(6,6,20)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,4,2))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,4,2))),lintype{II})
     else
        plot(axw,squeeze(squeeze(AmatNem(:,4,2))),lintype{II})
    end
       plot_properties

    ylabel([ylab,'_{42}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(6,6,21)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,4,3))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,4,3))),lintype{II})
     else
        plot(axw,squeeze(squeeze(AmatNem(:,4,3))),lintype{II})
    end
       plot_properties

    ylabel([ylab,'_{43}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(6,6,22)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,4,4))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,4,4))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,4,4))),lintype{II})
    end
    plot_properties

    ylabel([ylab,'_{44}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(6,6,25)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,5,1))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,5,1))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,5,1))),lintype{II})
    end
     plot_properties

    ylabel([ylab,'_{51}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(6,6,26)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,5,2))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,5,2))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,5,2))),lintype{II})
    end
     plot_properties

    ylabel([ylab,'_{52}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(6,6,27)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,5,3))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,5,3))),lintype{II})
     else
        plot(axw,squeeze(squeeze(AmatNem(:,5,3))),lintype{II})
    end
        plot_properties

    ylabel([ylab,'_{53}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(6,6,28)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,5,4))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,5,4))),lintype{II})
     else
        plot(axw,squeeze(squeeze(AmatNem(:,5,4))),lintype{II})
    end
     plot_properties

    ylabel([ylab,'_{54}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(6,6,29)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,5,5))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,5,5))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,5,5))),lintype{II})
    end
    ylabel([ylab,'_{55}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    plot_properties

    subplot(6,6,31)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,6,1))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,6,1))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,6,1))),lintype{II})
    end
    plot_properties

    ylabel([ylab,'_{61}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(6,6,32)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,6,2))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,6,2))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,6,2))),lintype{II})
    end
    plot_properties
    ylabel([ylab,'_{62}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(6,6,33)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,6,3))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,6,3))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,6,3))),lintype{II})
    end
    plot_properties

    ylabel([ylab,'_{63}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(6,6,34)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,6,4))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,6,4))),lintype{II})
     else
        plot(axw,squeeze(squeeze(AmatNem(:,6,4))),lintype{II})
    end
       plot_properties

    ylabel([ylab,'_{64}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(6,6,35)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,6,5))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,6,5))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,6,5))),lintype{II})
    end
     plot_properties

    ylabel([ylab,'_{65}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(6,6,36)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,6,6))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,6,6))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,6,6))),lintype{II})
    end
      plot_properties

    ylabel([ylab,'_{66}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
elseif ID_MODES_PLOT==2
    subplot(3,3,1)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,1,1))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,1,1))),lintype{II})
     else
        plot(axw,squeeze(squeeze(AmatNem(:,1,1))),lintype{II})
    end
         plot_properties

    ylabel([ylab,'_{11}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(3,3,4)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,3,3))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,3,3))),lintype{II})
      else
        plot(axw,squeeze(squeeze(AmatNem(:,3,3))),lintype{II})
    end
        plot_properties

    ylabel([ylab,'_{33}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(3,3,5)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,3,1))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,3,1))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,3,1))),lintype{II})
    end
    plot_properties

    ylabel([ylab,'_{31}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(3,3,7)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,5,5))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,5,5))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,5,5))),lintype{II})
    end
      plot_properties

    ylabel([ylab,'_{55}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(3,3,8)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,5,1))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,5,1))),lintype{II})
     else
        plot(axw,squeeze(squeeze(AmatNem(:,5,1))),lintype{II})
    end
     plot_properties

    ylabel([ylab,'_{51}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(3,3,9)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,5,3))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,5,3))),lintype{II})
     else
        plot(axw,squeeze(squeeze(AmatNem(:,5,3))),lintype{II})
    end
      plot_properties

    ylabel([ylab,'_{53}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
elseif ID_MODES_PLOT==3
    subplot(3,2,1)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,1,1))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,1,1))),lintype{II})
     else
        plot(axw,squeeze(squeeze(AmatNem(:,1,1))),lintype{II})
    end
        plot_properties

    ylabel([ylab,'_{11}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(3,2,2)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,2,2))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,2,2))),lintype{II})
     else
        plot(axw,squeeze(squeeze(AmatNem(:,2,2))),lintype{II})
    end
      plot_properties

    ylabel([ylab,'_{22}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(3,2,3)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,3,3))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,3,3))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,3,3))),lintype{II})
    end
      plot_properties

    ylabel([ylab,'_{33}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(3,2,4)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,4,4))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,4,4))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,4,4))),lintype{II})
    end
      plot_properties

    ylabel([ylab,'_{44}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(3,2,5)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,5,5))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,5,5))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,5,5))),lintype{II})
    end
     plot_properties

    ylabel([ylab,'_{55}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(3,2,6)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(AMatH(:,6,6))),HySmb,...
            axw,squeeze(squeeze(AmatNem(:,6,6))),lintype{II})
    else
        plot(axw,squeeze(squeeze(AmatNem(:,6,6))),lintype{II})
    end
      plot_properties

    ylabel([ylab,'_{66}']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
end
end

function fun_plot_excitationForcesRAO(axwH,axw,FdatH,FdatNem,HySmb,II,...
    lintype,Idlimitaxes,xlimval,xlab,IDHYDROSTAR,IdbetaPlotH,IdbetaPlot,ylab,ylabN,ID_MODES_PLOT)
if ID_MODES_PLOT==1
subplot(6,2,1)
if II==1 && IDHYDROSTAR==1
    plot(axwH,squeeze(squeeze(FdatH(:,1,IdbetaPlotH))),HySmb,...
        axw,squeeze(squeeze(FdatNem(:,1,IdbetaPlot))),lintype{II})
else
    plot(axw,squeeze(squeeze(FdatNem(:,1,IdbetaPlot))),lintype{II})
end
 plot_properties

ylabel([ylab,'_{1}',ylabN]);
hold on;
if Idlimitaxes==1
    xlim([xlimval]);
end
xlabel(xlab)
subplot(6,2,2)
if II==1 && IDHYDROSTAR==1
    plot(axwH,squeeze(squeeze(FdatH(:,7,IdbetaPlotH))),HySmb,...
        axw,squeeze(squeeze(FdatNem(:,7,IdbetaPlot))),lintype{II})
else
    plot(axw,squeeze(squeeze(FdatNem(:,7,IdbetaPlot))),lintype{II})
end
  plot_properties

ylabel(['arg(',ylab(2:end),'_{1}',ylabN(1:end-1),')']);
hold on;
if Idlimitaxes==1
    xlim([xlimval]);
end
xlabel(xlab)

subplot(6,2,3)
if II==1 && IDHYDROSTAR==1
    plot(axwH,squeeze(squeeze(FdatH(:,2,IdbetaPlotH))),HySmb,...
        axw,squeeze(squeeze(FdatNem(:,2,IdbetaPlot))),lintype{II})
else
    plot(axw,squeeze(squeeze(FdatNem(:,2,IdbetaPlot))),lintype{II})
end
 plot_properties

ylabel([ylab,'_{2}',ylabN]);
hold on;
if Idlimitaxes==1
    xlim([xlimval]);
end
xlabel(xlab)
subplot(6,2,4)
if II==1 && IDHYDROSTAR==1
    plot(axwH,squeeze(squeeze(FdatH(:,8,IdbetaPlotH))),HySmb,...
        axw,squeeze(squeeze(FdatNem(:,8,IdbetaPlot))),lintype{II})
else
    plot(axw,squeeze(squeeze(FdatNem(:,8,IdbetaPlot))),lintype{II})
end
  plot_properties

ylabel(['arg(',ylab(2:end),'_{1}',ylabN(1:end-1),')']);
hold on;
if Idlimitaxes==1
    xlim([xlimval]);
end
xlabel(xlab)

subplot(6,2,5)
if II==1 && IDHYDROSTAR==1
    plot(axwH,squeeze(squeeze(FdatH(:,3,IdbetaPlotH))),HySmb,...
        axw,squeeze(squeeze(FdatNem(:,3,IdbetaPlot))),lintype{II})
else
    plot(axw,squeeze(squeeze(FdatNem(:,3,IdbetaPlot))),lintype{II})
end
 plot_properties

ylabel([ylab,'_{3}',ylabN]);
hold on;
if Idlimitaxes==1
    xlim([xlimval]);
end
xlabel(xlab)
subplot(6,2,6)
if II==1 && IDHYDROSTAR==1
    plot(axwH,squeeze(squeeze(FdatH(:,9,IdbetaPlotH))),HySmb,...
        axw,squeeze(squeeze(FdatNem(:,9,IdbetaPlot))),lintype{II})
else
    plot(axw,squeeze(squeeze(FdatNem(:,9,IdbetaPlot))),lintype{II})
end
  plot_properties

ylabel(['arg(',ylab(2:end),'_{3}',ylabN(1:end-1),')']);
hold on;
if Idlimitaxes==1
    xlim([xlimval]);
end
xlabel(xlab)

subplot(6,2,7)
if II==1 && IDHYDROSTAR==1
    plot(axwH,squeeze(squeeze(FdatH(:,4,IdbetaPlotH))),HySmb,...
        axw,squeeze(squeeze(FdatNem(:,4,IdbetaPlot))),lintype{II})
else
    plot(axw,squeeze(squeeze(FdatNem(:,4,IdbetaPlot))),lintype{II})
end
 plot_properties

ylabel([ylab,'_{4}',ylabN]);
hold on;
if Idlimitaxes==1
    xlim([xlimval]);
end
xlabel(xlab)
subplot(6,2,8)
if II==1 && IDHYDROSTAR==1
    plot(axwH,squeeze(squeeze(FdatH(:,10,IdbetaPlotH))),HySmb,...
        axw,squeeze(squeeze(FdatNem(:,10,IdbetaPlot))),lintype{II})
else
    plot(axw,squeeze(squeeze(FdatNem(:,10,IdbetaPlot))),lintype{II})
end
   plot_properties

ylabel(['arg(',ylab(2:end),'_{4}',ylabN(1:end-1),')']);
hold on;
if Idlimitaxes==1
    xlim([xlimval]);
end
xlabel(xlab)

subplot(6,2,9)
if II==1 && IDHYDROSTAR==1
    plot(axwH,squeeze(squeeze(FdatH(:,5,IdbetaPlotH))),HySmb,...
        axw,squeeze(squeeze(FdatNem(:,5,IdbetaPlot))),lintype{II})
else
    plot(axw,squeeze(squeeze(FdatNem(:,5,IdbetaPlot))),lintype{II})
end
 plot_properties

ylabel([ylab,'_{5}',ylabN]);
hold on;
if Idlimitaxes==1
    xlim([xlimval]);
end
xlabel(xlab)
subplot(6,2,10)
if II==1 && IDHYDROSTAR==1
    plot(axwH,squeeze(squeeze(FdatH(:,11,IdbetaPlotH))),HySmb,...
        axw,squeeze(squeeze(FdatNem(:,11,IdbetaPlot))),lintype{II})
else
    plot(axw,squeeze(squeeze(FdatNem(:,11,IdbetaPlot))),lintype{II})
end
 plot_properties

ylabel(['arg(',ylab(2:end),'_{5}',ylabN(1:end-1),')']);
hold on;
if Idlimitaxes==1
    xlim([xlimval]);
end
xlabel(xlab)

subplot(6,2,11)
if II==1 && IDHYDROSTAR==1
    plot(axwH,squeeze(squeeze(FdatH(:,6,IdbetaPlotH))),HySmb,...
        axw,squeeze(squeeze(FdatNem(:,6,IdbetaPlot))),lintype{II})
else
    plot(axw,squeeze(squeeze(FdatNem(:,6,IdbetaPlot))),lintype{II})
end
 plot_properties

ylabel([ylab,'_{6}',ylabN]);
hold on;
if Idlimitaxes==1
    xlim([xlimval]);
end
xlabel(xlab)
subplot(6,2,12)
if II==1 && IDHYDROSTAR==1
    plot(axwH,squeeze(squeeze(FdatH(:,12,IdbetaPlotH))),HySmb,...
        axw,squeeze(squeeze(FdatNem(:,12,IdbetaPlot))),lintype{II})
else
    plot(axw,squeeze(squeeze(FdatNem(:,12,IdbetaPlot))),lintype{II})
end
  plot_properties

ylabel(['arg(',ylab(2:end),'_{6}',ylabN(1:end-1),')']);
hold on;
if Idlimitaxes==1
    xlim([xlimval]);
end
xlabel(xlab)

elseif ID_MODES_PLOT==2
    subplot(3,2,1)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(FdatH(:,1,IdbetaPlotH))),HySmb,...
            axw,squeeze(squeeze(FdatNem(:,1,IdbetaPlot))),lintype{II})
       else
        plot(axw,squeeze(squeeze(FdatNem(:,1,IdbetaPlot))),lintype{II})
    end
      plot_properties

    ylabel([ylab,'_{1}',ylabN]);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(3,2,2)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(FdatH(:,7,IdbetaPlotH))),HySmb,...
            axw,squeeze(squeeze(FdatNem(:,7,IdbetaPlot))),lintype{II})
    else
        plot(axw,squeeze(squeeze(FdatNem(:,7,IdbetaPlot))),lintype{II})
    end
       plot_properties
   ylabel(['arg(',ylab(2:end),'_{1}',ylabN(1:end-1),')']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(3,2,3)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(FdatH(:,3,IdbetaPlotH))),HySmb,...
            axw,squeeze(squeeze(FdatNem(:,3,IdbetaPlot))),lintype{II})
    else
        plot(axw,squeeze(squeeze(FdatNem(:,3,IdbetaPlot))),lintype{II})
    end
      plot_properties;
    ylabel([ylab,'_{3}',ylabN]);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(3,2,4)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(FdatH(:,9,IdbetaPlotH))),HySmb,...
            axw,squeeze(squeeze(FdatNem(:,9,IdbetaPlot))),lintype{II})
    else
        plot(axw,squeeze(squeeze(FdatNem(:,9,IdbetaPlot))),lintype{II})
    end
      plot_properties
   ylabel(['arg(',ylab(2:end),'_{3}',ylabN(1:end-1),')']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

     subplot(3,2,5)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(FdatH(:,5,IdbetaPlotH))),HySmb,...
            axw,squeeze(squeeze(FdatNem(:,5,IdbetaPlot))),lintype{II})
    else
        plot(axw,squeeze(squeeze(FdatNem(:,5,IdbetaPlot))),lintype{II})
    end
    plot_properties
    ylabel([ylab,'_{5}',ylabN]);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(3,2,6)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(FdatH(:,11,IdbetaPlotH))),HySmb,...
            axw,squeeze(squeeze(FdatNem(:,11,IdbetaPlot))),lintype{II})
    else
        plot(axw,squeeze(squeeze(FdatNem(:,11,IdbetaPlot))),lintype{II})
    end
    plot_properties
   ylabel(['arg(',ylab(2:end),'_{5}',ylabN(1:end-1),')']);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    elseif ID_MODES_PLOT==3
    subplot(3,2,1)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(FdatH(:,1,IdbetaPlotH))),HySmb,...
            axw,squeeze(squeeze(FdatNem(:,1,IdbetaPlot))),lintype{II})
    else
        plot(axw,squeeze(squeeze(FdatNem(:,1,IdbetaPlot))),lintype{II})
    end
    plot_properties
    ylabel([ylab,'_{1}',ylabN]);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(3,2,2)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(FdatH(:,2,IdbetaPlotH))),HySmb,...
            axw,squeeze(squeeze(FdatNem(:,2,IdbetaPlot))),lintype{II})
    else
        plot(axw,squeeze(squeeze(FdatNem(:,7,IdbetaPlot))),lintype{II})
    end
    plot_properties
   ylabel([ylab,'_{2}',ylabN]);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

    subplot(3,2,3)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(FdatH(:,3,IdbetaPlotH))),HySmb,...
            axw,squeeze(squeeze(FdatNem(:,3,IdbetaPlot))),lintype{II})
    else
        plot(axw,squeeze(squeeze(FdatNem(:,3,IdbetaPlot))),lintype{II})
    end
    plot_properties
    ylabel([ylab,'_{3}',ylabN]);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(3,2,4)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(FdatH(:,4,IdbetaPlotH))),HySmb,...
            axw,squeeze(squeeze(FdatNem(:,4,IdbetaPlot))),lintype{II})
    else
        plot(axw,squeeze(squeeze(FdatNem(:,4,IdbetaPlot))),lintype{II})
    end
         plot_properties

    ylabel([ylab,'_{4}',ylabN]);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

     subplot(3,2,5)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(FdatH(:,5,IdbetaPlotH))),HySmb,...
            axw,squeeze(squeeze(FdatNem(:,5,IdbetaPlot))),lintype{II})
    else
        plot(axw,squeeze(squeeze(FdatNem(:,5,IdbetaPlot))),lintype{II})
    end
     plot_properties

    ylabel([ylab,'_{5}',ylabN]);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)
    subplot(3,2,6)
    if II==1 && IDHYDROSTAR==1
        plot(axwH,squeeze(squeeze(FdatH(:,6,IdbetaPlotH))),HySmb,...
            axw,squeeze(squeeze(FdatNem(:,6,IdbetaPlot))),lintype{II})
    else
        plot(axw,squeeze(squeeze(FdatNem(:,6,IdbetaPlot))),lintype{II})
    end
        plot_properties

   ylabel([ylab,'_{6}',ylabN]);
    hold on;
    if Idlimitaxes==1
        xlim([xlimval]);
    end
    xlabel(xlab)

end
end




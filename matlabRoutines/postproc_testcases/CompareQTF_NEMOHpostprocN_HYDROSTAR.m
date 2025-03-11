%code for comparing the QTF results between NEMOH & Hydrostar
%data parameters are set in data_param;
%make sure the testcases data available after run.
%R. Kurnia, 2022 ECN

close all;
clc;clear;

[pathstr,~,~] = fileparts(mfilename('fullpath'));
cd (pathstr);
addpath(genpath(pathstr))
%%%
TestCaseNo=1;
if TestCaseNo==1        %8b_QTF_Cylinder
    dataparamQTF_1;
    qtftype='M';
    betaID=[0 0];
    maxValQTF=30;
    Idwlim=1; wlimval=[0,pi];
elseif TestCaseNo==2    %9_QTF_OC4_Semisubmersible
    dataparamQTF_2;
    qtftype='M';
    betaID=[0 30];      % choose beta1, beta2 !
     maxValQTF=20;
     Idwlim=1; wlimval=[0,pi];
elseif TestCaseNo==3     %10_QTF_SOFTWIND
    dataparamQTF_3;
    qtftype='M';
    betaID=[0 30];       % choose beta1, beta2 !
     maxValQTF=20;%20;5;1000
     Idwlim=1; wlimval=[0,pi];
elseif TestCaseNo==4     %10_QTF_SOFTWIND_FS
    dataparamQTF_4;
    qtftype='P';
    betaID=[0 0];        % choose beta1, beta2 !
    maxValQTF=40;%40;10;2000
    Idwlim=1; wlimval=[0,2];
elseif TestCaseNo==5     %11_QTF_OC3_Hywind
    dataparamQTF_5;
    qtftype='M';
    betaID=[0 0];        % choose beta1, beta2 !
    maxValQTF=6;%6,1,50
    Idwlim=1; wlimval=[0,pi];
end
 DOF=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Idwvar=0;%0=w rad/s,1=f 1/s,2=T s;
%Idwlim=0; wlimval=[0.06,0.2];
%Idwlim=1; wlimval=[0.05,1.5];
%Idwlim=1; wlimval=[0.05,1.5];
Idaxisij=0;
IDHYDROSTAR=1;
NHEADERHS=22;%number of header line in the Hydrostar data
betaIDHS=circshift(betaID,1);
DOFdatNem=[1 1 1 1 1 1];% Surge, Heave, pitch
shiftdw1=-1;% mean drift
shiftdw2=-3;%surge/pitch natural freq
IdCalTerms=0; %0: plot total QTF data from NEMOH post_processing output
%1: plot total QTF from all the saved QTF terms
%1: for debugging then choose the terms to be plotted below
IdNemDuokTerm=[1 1 1 1 1 1]*1; %1 computed 0 not,
%[DUOK_QUADVEL ...] look terms in NEMOH
% all zero mean duok terms are not computed
IdNemHasboTerm=[1 1 1 1 1 0]*1; %1 computed 0 not,
IdNemHasfsTerm=[1 1 1 1 1 0 0 0 0 ]*1; %1 computed 0 not,
%IdNemHasfsTerm=[1 1 1 0 0 1 1 1 1]; %1 computed 0 not,
% Index 1-3: QFD1
% Index 4-5: QFD2 direct formulation dzz=k^2
% Index 6-9: QFD2 Green Formulation
IdNemAsympTerm=[1 1]*1; %1 computed 0 not,
% Index 1-2: QFD_ASYMP
ShowLegend=0;

IdNormDiff =0; %1 normalized by max value, 0 only difference H-N
IdContval=1;
IdContvaldif=1;
IdMinMaxQTF=1;
MinMaxQTF=[0 maxValQTF];%QTF/rho/g
IdMinMaxQTFdiff=1;
MinMaxQTFdiff=[-max(MinMaxQTF),max(MinMaxQTF)]/10;

IdMinMaxQTFReIm=1;
MinMaxReIm=[-max(MinMaxQTF),max(MinMaxQTF)];
MinMaxAbs=[0,MinMaxReIm(2)];


if DOF==1
  %  Contval=[10,20,40,60,100];
   % Contval=[10,20,40,60];
   % Contval=[10,20,30];
      Contval=[2,5,10];
   %Contval=[0.2:0.1:1,1.5:0.5:5];
   %Contval=[0.2,0.5,1,2:1:10,12:2:26];
  % Contval=[0.2:0.1:0.5,0.6:0.2:1,1.5:0.5:5];
  % Contval=[10,20];
   %Contval=[2,5,10,15];
elseif DOF==2
    %Contval=[5,10,20,40,50];
    %Contval=[2,5,10,20];
    Contval=[2,5,8];

elseif DOF==3
  %  Contval=[5,10,20,40];
  %   Contval=[5,10,15];
  % Contval=[2,5,10,20];
    Contval=[1,2,4];
   % Contval=[0.5,1,2];
  %  Contval=[0.5,1,1.5];
elseif (DOF==4)
    Contval=[100,200,300];
elseif (DOF==5)
  %   Contval=[500,1000,1500];
    Contval=[100,200,400,800,1000];
   %  Contval=[100,200,400,600];
  %  Contval=[10,20,40,80];
   % Contval=[10,20,40];
elseif (DOF==6)
    Contval=[200,500,1000];
end

Contvaldif=floor(Contval/10);%[20,200,500];
Contvaldif(Contvaldif==0)=1;

%Contvaldif=[20,50,100,500];
% Contvaldif=[2,5,10,30,60,100];
% Contvaldif=[10,20,50,100,200];

rho=1025;
g=9.81;

if IdCalTerms==0
    [QdatTotR,QdatTotI,QdatTotMod,w1]=fun_readNEMOHQTFOutput_N...
        (['./',NEMOH_projdir,'/results/QTF/OUT_QTF',qtftype,'_N.dat'],...
        NWdatNEM,DOFdatNem,DOF,betaID,NbetaData,qtftype,SwitchBiDir); %data output is already normalized by rho g
else
    NDterm=length(IdNemDuokTerm);
    QdatTotR=0;QdatTotI=0;
    for i=1:NDterm
        if (IdNemDuokTerm(i)==1)
            [QdatTotR_I,QdatTotI_I,~,w1]=fun_readNEMOHQTFOutputTerms_N...
                (['./',NEMOH_projdir,'/results/QTF/QTF',qtftype,'_DUOK_term_',num2str(i),'.dat'],...
                NWdatNEM,DOFdatNem,DOF,betaID,NbetaData,qtftype,rho,g); %data output is already normalized by rho g
            QdatTotR=QdatTotR+QdatTotR_I;QdatTotI=QdatTotI+QdatTotI_I;
        end
    end
    NHterm=length(IdNemHasboTerm);
    for i=1:NHterm
        if (IdNemHasboTerm(i)==1)
            [QdatTotR_I,QdatTotI_I,~,w1]=fun_readNEMOHQTFOutputTerms_N...
                (['./',NEMOH_projdir,'/results/QTF/QTF',qtftype,'_HASBO_term_',num2str(i),'.dat'],...
                NWdatNEM,DOFdatNem,DOF,betaID,NbetaData,qtftype,rho,g); %data output is already normalized by rho g
            QdatTotR=QdatTotR+QdatTotR_I;QdatTotI=QdatTotI+QdatTotI_I;
        end
    end

    NHterm=length(IdNemHasfsTerm);
    for i=1:NHterm
        if (IdNemHasfsTerm(i)==1)
            [QdatTotR_I,QdatTotI_I,~,w1]=fun_readNEMOHQTFOutputTerms_N...
                (['./',NEMOH_projdir,'/results/QTF/QTF',qtftype,'_HASFS_term_',num2str(i),'.dat'],...
                NWdatNEM,DOFdatNem,DOF,betaID,NbetaData,qtftype,rho,g); %data output is already normalized by rho g
            QdatTotR=QdatTotR+QdatTotR_I;QdatTotI=QdatTotI+QdatTotI_I;
        end
    end

    NHterm=length(IdNemAsympTerm);
    for i=1:NHterm
        if (IdNemAsympTerm(i)==1)
            [QdatTotR_I,QdatTotI_I,~,w1]=fun_readNEMOHQTFOutputTerms_N...
                (['./',NEMOH_projdir,'/results/QTF/QTF',qtftype,'_ASYMP_term_',num2str(i),'.dat'],...
                NWdatNEM,DOFdatNem,DOF,betaID,NbetaData,qtftype,rho,g); %data output is already normalized by rho g
            QdatTotR=QdatTotR+QdatTotR_I;QdatTotI=QdatTotI+QdatTotI_I;
        end
    end

    QdatTotMod=sqrt(QdatTotR.^2+QdatTotI.^2);
end

w1lab='\omega_1 [rad/s]';
w2lab='\omega_2 [rad/s]';

if Idwvar==1
    w1=w1./2/pi;
    w1lab='f_1 [Hz]';
    w2lab='f_2 [Hz]';
elseif Idwvar==2
    w1=2*pi./w1;
    w1lab='T_1 [s]';
    w2lab='T_2 [s]';
end

if Idwlim==1
    Idwi=closest(w1,wlimval(1));
    Idwf=closest(w1,wlimval(2));
else
    Idwi=1;
    Idwf=length(w1);
end
w1=w1(Idwi:Idwf);
QdatTotR=QdatTotR(Idwi:Idwf,Idwi:Idwf);
QdatTotI=QdatTotI(Idwi:Idwf,Idwi:Idwf);
QdatTotMod=QdatTotMod(Idwi:Idwf,Idwi:Idwf);

[ww1,ww2]=meshgrid(w1,w1);


if DOF==1
    Qtfdof='Fx';
elseif DOF==2
    Qtfdof='Fy';
elseif DOF==3
    Qtfdof='Fz';
elseif DOF==4
    Qtfdof='Mx';
elseif DOF==5
    Qtfdof='My';
elseif DOF==6
    Qtfdof='Mz';
end

if IDHYDROSTAR==1
    if NbetaDataHydrodatar>1
        if strcmpi(qtftype,'M')
                [HSQTFTotRe,HSQTFTotIm,wH1]=fun_readHydrostarMQTFOutput(['./',hydrostarDir,'/MQTF',Qtfdof,'.rao'],NHEADERHS,NWHS,betaIDHS,NbetaDataHydrodatar,qtftype);
        elseif strcmpi(qtftype,'P')
                [HSQTFTotRe,HSQTFTotIm,wH1]=fun_readHydrostarMQTFOutput(['./',hydrostarDir,'/MSPG',Qtfdof,'.rao'],NHEADERHS,NWHS,betaIDHS,NbetaDataHydrodatar,qtftype);
        end
    else
        if strcmpi(qtftype,'M')
        [HSQTFTotRe,HSQTFTotIm,wH1]=fun_readHydrostarQTFOutput(['./',hydrostarDir,'/QTF',Qtfdof,'.rao'],NHEADERHS,NWHS,qtftype);
        elseif strcmpi(qtftype,'P')
        [HSQTFTotRe,HSQTFTotIm,wH1]=fun_readHydrostarQTFOutput(['./',hydrostarDir,'/SPG',Qtfdof,'.rao'],NHEADERHS,NWHS,qtftype);
        end
    end
    if Idwvar==1
        wH1=wH1./2/pi;
    elseif Idwvar==2
        wH1=2*pi./wH1;
    end

    if Idwlim==1
        IdwHi=closest(wH1,wlimval(1));
        IdwHf=closest(wH1,wlimval(2));
    else
        IdwHi=1;
        IdwHf=length(wH1);
    end
    wH1=wH1(IdwHi:IdwHf);
    [wwH1,wwH2]=meshgrid(wH1,wH1);

    if length(w1)==length(wH1)
    HSQTFTotRe=HSQTFTotRe(IdwHi:IdwHf,IdwHi:IdwHf)./rho/g;
    HSQTFTotIm=HSQTFTotIm(IdwHi:IdwHf,IdwHi:IdwHf)./rho/g;
    else
    HSQTFTotRe=interp2(wwH1,wwH2,HSQTFTotRe(IdwHi:IdwHf,IdwHi:IdwHf),ww1,ww2)./rho/g;
    HSQTFTotIm=interp2(wwH1,wwH2,HSQTFTotIm(IdwHi:IdwHf,IdwHi:IdwHf),ww1,ww2)./rho/g;
    wwH1=ww1;
    wwH2=ww2;
    end

    HSQTFTotMod=sqrt(HSQTFTotIm.^2+HSQTFTotRe.^2);


    if  IdMinMaxQTFReIm==0
    minTotR=min(min(min(QdatTotR)),min(min(HSQTFTotRe)));
    maxTotR=max(max(max(QdatTotR)),max(max(HSQTFTotRe)));
    minTotI=min(min(min(QdatTotI)),min(min(HSQTFTotIm)));
    maxTotI=max(max(max(QdatTotI)),max(max(HSQTFTotIm)));
    minTotMod=min(min(min(QdatTotMod)),min(min(HSQTFTotMod)));
    maxTotMod=max(max(max(QdatTotMod)),max(max(HSQTFTotMod)));
    else
    minTotI=MinMaxReIm(1);maxTotI=MinMaxReIm(2);
    minTotR=minTotI(1);maxTotR=maxTotI;
    minTotMod=MinMaxAbs(1);maxTotMod=MinMaxAbs(2);
    end

    NormalizedAbsDifTOTR=abs((QdatTotR-HSQTFTotRe))./max(max(abs(HSQTFTotRe)));
    NormalizedAbsDifTOTR(isnan( NormalizedAbsDifTOTR))=0;
    NormalizedAbsDifTOTI=abs((QdatTotI-HSQTFTotIm))./max(max(abs(HSQTFTotIm)));
    NormalizedAbsDifTOTI(isnan( NormalizedAbsDifTOTI))=0;
    NormalizedAbsDifTOTMOD=abs((QdatTotMod-HSQTFTotMod))./max(max(abs(HSQTFTotMod)));
    NormalizedAbsDifTOTMOD(isnan( NormalizedAbsDifTOTMOD))=0;

    QdiagR1_HS=diag(HSQTFTotRe,shiftdw1);
    QdiagI1_HS=diag(HSQTFTotIm,shiftdw1);
    QdiagA1_HS=diag(sqrt(HSQTFTotRe.^2+HSQTFTotIm.^2),shiftdw1);
    QdiagP1_HS=diag(atan(HSQTFTotIm./HSQTFTotRe),shiftdw1)./pi;

    QdiagR2_HS=diag(HSQTFTotRe,shiftdw2);
    QdiagI2_HS=diag(HSQTFTotIm,shiftdw2);
    QdiagA2_HS=diag(sqrt(HSQTFTotRe.^2+HSQTFTotIm.^2),shiftdw2);
    QdiagP2_HS=diag(atan(HSQTFTotIm./HSQTFTotRe),shiftdw2)./pi;
end

QdiagR1=diag(QdatTotR,shiftdw1);
QdiagI1=diag(QdatTotI,shiftdw1);
QdiagA1=diag(sqrt(QdatTotR.^2+QdatTotI.^2),shiftdw1);
QdiagP1=diag(atan(QdatTotI./QdatTotR),shiftdw1)./pi;
QdiagP1(QdiagP1>0.3)=-QdiagP1(QdiagP1>0.3);

QdiagR2=diag(QdatTotR,shiftdw2);
QdiagI2=diag(QdatTotI,shiftdw2);
QdiagA2=diag(sqrt(QdatTotR.^2+QdatTotI.^2),shiftdw2);
QdiagP2=diag(atan(QdatTotI./QdatTotR),shiftdw2)./pi;
QdiagP2(QdiagP2>0.3)=-QdiagP2(QdiagP2>0.3);

freqnorm1=w1(abs(shiftdw1)+1:end);
freqnorm2=w1(2:end-(abs(shiftdw2)-1));
dww=w1(2)-w1(1);
%% line plot QTF diag
if strcmpi(qtftype,'P')
    qtflabel='QTF+';
else
    qtflabel='QTF-';
end
figure;
subplot(2,2,1)
if IDHYDROSTAR==1
    plot(freqnorm1,QdiagA1_HS,'b',...
        freqnorm1,QdiagA1,'--r');
    legend('HYDROSTAR','NEMOH')
else
    plot(freqnorm1,QdiagA1,'--r');
    legend('NEMOH')
end
if Idwvar==0
    title(['Diag. ', qtflabel,' for \Deltaw=',num2str(abs(shiftdw1)*dww,3)])
elseif Idwvar==1
    title(['Diag. ', qtflabel,' for \Deltaf=',num2str(abs(shiftdw1)*dww,3)])
else
    title(['Diag. ', qtflabel,' for \DeltaT=',num2str(abs(shiftdw1)*dww,3)])
end
%ylim([0 0.8]);xlim([0 max(freqnorm1)])
xlabel(w1lab);
ylabel(['|', qtflabel,'|/\rhog'])
plot_properties;

subplot(2,2,2)

if IDHYDROSTAR==1
    plot(freqnorm1,-QdiagP1_HS,'b',...
        freqnorm1,-wrapToPi(QdiagP1),'--r');
    legend('HYDROSTAR','NEMOH')
else
    plot(freqnorm1,-wrapToPi(QdiagP1),'--r');
    legend('NEMOH')
end

if Idwvar==0
    title([qtflabel,' for \Deltaw=',num2str(abs(shiftdw1)*dww,3)])
elseif Idwvar==1
    title([qtflabel,' for \Deltaf=',num2str(abs(shiftdw1)*dww,3)])
else
    title([qtflabel,' for \DeltaT=',num2str(abs(shiftdw1)*dww,3)])
end
%ylim([0 0.8]);xlim([0 max(freqnorm1)])
xlabel(w1lab);
ylabel(['-\angle ', qtflabel,'/\pi'])
plot_properties;

subplot(2,2,3)

if IDHYDROSTAR==1
    plot(freqnorm2,QdiagA2_HS,'b',...
        freqnorm2,QdiagA2,'--r');
    legend('HYDROSTAR','NEMOH')
else
    plot(freqnorm2,QdiagA2,'--r');
    legend('NEMOH')
end

if Idwvar==0
    title(['Diag. ', qtflabel,' for \Deltaw=',num2str(abs(shiftdw2)*dww,3)])
elseif Idwvar==1
    title(['Diag. ', qtflabel,' for \Deltaf=',num2str(abs(shiftdw2)*dww,3)])
else
    title(['Diag. ', qtflabel,' for \DeltaT=',num2str(abs(shiftdw2)*dww,3)])
end
%ylim([0 0.8]);xlim([0 max(freqnorm1)])
xlabel(w1lab);
ylabel(['|', qtflabel,'|/\rhog'])
plot_properties;

subplot(2,2,4)

if IDHYDROSTAR==1
    plot(freqnorm2,-QdiagP2_HS,'b',...
        freqnorm2,-wrapToPi(QdiagP2),'--r');
    legend('HYDROSTAR','NEMOH')
else
    plot(freqnorm2,-wrapToPi(QdiagP2),'--r');
    legend('NEMOH')
end

if Idwvar==0
    title([qtflabel,' for \Deltaw=',num2str(abs(shiftdw2)*dww,3)])
elseif Idwvar==1
    title([qtflabel,' for \Deltaf=',num2str(abs(shiftdw2)*dww,3)])
else
    title([qtflabel,' for \DeltaT=',num2str(abs(shiftdw2)*dww,3)])
end
%ylim([0 0.8]);xlim([0 max(freqnorm1)])
xlabel(w1lab);
ylabel(['-\angle ', qtflabel,'/\pi'])
plot_properties;


figure;
subplot(2,1,1)
if IDHYDROSTAR==1
    plot(freqnorm1,QdiagR1_HS,'-b',freqnorm1,QdiagI1_HS,'-.b',freqnorm1,QdiagR1,'--r',freqnorm1,QdiagI1,':r');
    if ShowLegend==1
        legend('HYDROSTAR Real','-HYDROSTAR Imag.','NEMOH Real','NEMOH Imag.')
    end
else
    plot(freqnorm1,QdiagR1,'--r',freqnorm1,QdiagI1,':r');
    if ShowLegend==1
        legend('NEMOH Real','NEMOH Imag.')
    end
end
if Idwvar==0
    xlabel('$\omega [rad/s] $','interpreter','latex');
elseif Idwvar==1
    xlabel('$f[Hz] $','interpreter','latex');
else
    xlabel('$T [s] $','interpreter','latex');
end

if Idwvar==0
    title([qtflabel,' for \Deltaw=',num2str(abs(shiftdw1)*dww,3)])
elseif Idwvar==1
    title([qtflabel,' for \Deltaf=',num2str(abs(shiftdw1)*dww,3)])
else
    title([qtflabel,' for \DeltaT=',num2str(abs(shiftdw1)*dww,3)])
end

%ylim([0 0.8]);
xlim([min(freqnorm1) max(freqnorm1)])
ylabel([qtflabel,'/\rhog'])
plot_properties;
subplot(2,1,2)
if IDHYDROSTAR==1
    plot(freqnorm2,QdiagR2_HS,'-b',freqnorm2,QdiagI2_HS,'-.b',freqnorm2,QdiagR2,'--r',freqnorm2,QdiagI2,':r');

else
    plot(freqnorm2,QdiagR2,'--r',freqnorm2,QdiagI2,':r');

end
if Idwvar==0
    title([qtflabel,' for \Deltaw=',num2str(abs(shiftdw2)*dww,3)])
elseif Idwvar==1
    title([qtflabel,' for \Deltaf=',num2str(abs(shiftdw2)*dww,3)])
else
    title([qtflabel,' for \DeltaT=',num2str(abs(shiftdw2)*dww,3)])
end

if Idwvar==0
    xlabel('$\omega [rad/s] $','interpreter','latex');
elseif Idwvar==1
    xlabel('$f[Hz] $','interpreter','latex');
else
    xlabel('$T [s] $','interpreter','latex');
end
%ylim([0 0.8]);
xlim([min(freqnorm2) max(freqnorm2)])
ylabel([qtflabel,'/\rhog'])
plot_properties;


if IDHYDROSTAR==1
    %% [c,h]=contour plot QTF total
    if  IdMinMaxQTF==1
        minQMod=MinMaxQTF(1);
        maxQMod=MinMaxQTF(2);
    else
        minQMod=min(min(min(QdatTotMod)),min(min(HSQTFTotMod)));
        maxQMod=max(max(max(QdatTotMod)),max(max(HSQTFTotMod)));
    end
    [WW1,WW2]=meshgrid(ww1,ww2);
    figure

    diffQdatTotMod=HSQTFTotMod-QdatTotMod;
    subplot(1,3,2)
    surf(ww1,ww2,QdatTotMod,'edgecolor','none');
    ax=gca;
    hold(ax,'on');
    if IdContval==1
         funP_contour(ax,ww1,ww2,QdatTotMod,Contval,max(max(max(QdatTotMod))),'w',1.5);
    else
        funP_contour(ax,ww1,ww2,QdatTotMod,[],max(max(max(QdatTotMod))),'w',1.5);
    end
    hold(ax,'off');
    title(['NEMOH'])
    xlabel(w1lab);ylabel(w2lab);
    colormap('jet')
    caxis([minQMod, maxQMod])
     cb=colorbar;
     ylabel(cb,[qtflabel,'/\rhog'])
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    plot_properties;
    subplot(1,3,1)
    surf(ww1,ww2,HSQTFTotMod,'edgecolor','none');
    ax=gca;
    hold(ax,'on');
    if IdContval==1
        funP_contour(ax,ww1,ww2,HSQTFTotMod,Contval,max(max(max(HSQTFTotMod))),'w',1.5);
    else
        funP_contour(ax,ww1,ww2,HSQTFTotMod,[],max(max(max(HSQTFTotMod))),'w',1.5);
    end

    xlabel(w1lab);ylabel(w2lab);
    title(['HYDROSTAR'])
    caxis([minQMod, maxQMod])
    xlabel(w1lab);ylabel(w2lab);
    colormap('jet')
     cb=colorbar;
      ylabel(cb,[qtflabel,'/\rhog'])
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    plot_properties;
    subplot(1,3,3)
    if IdNormDiff==0
        surf(ww1,ww2,diffQdatTotMod,'edgecolor','none');
        ax=gca;
        hold(ax,'on');
        if IdContval==1
            funP_contour(ax,ww1,ww2,diffQdatTotMod,Contvaldif,max(max(max(diffQdatTotMod))),'w',1.5);
        else
            funP_contour(ax,ww1,ww2,diffQdatTotMod,[],max(max(max(diffQdatTotMod))),'w',1.5);
        end
    else
        surf(ww1,ww2,NormalizedAbsDifTOTMOD,'edgecolor','none');
        ax=gca;
        hold(ax,'on');
        if IdContval==1
            funP_contour(ax,ww1,ww2,NormalizedAbsDifTOTMOD,Contvaldif,max(max(max(NormalizedAbsDifTOTMOD))),'w',1.5);
        else
            funP_contour(ax,ww1,ww2,NormalizedAbsDifTOTMOD,[],max(max(max(NormalizedAbsDifTOTMOD))),'w',1.5);
        end
    end
     xlabel(w1lab);ylabel(w2lab);
    title(['Normalized. Abs. Diff.'])
    if IdNormDiff==0
        title(['Difference'])
    else
        title(['Normalized. Abs. Diff.'])
    end
    xlabel(w1lab);ylabel(w2lab);
    if  IdMinMaxQTFdiff==1
    caxis([MinMaxQTFdiff(1), MinMaxQTFdiff(2)])
    end
    colormap('jet');
      cb=colorbar;
   %  ylabel(cb,[qtflabel,'/\rhog'])
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    plot_properties;


    %% surf plot QTF total
    figure;
    subplot(3,3,1)
    surf(ww1,ww2,QdatTotR,'edgecolor','none')
    title(['NEMOH : Re Tot'])
    zlim([minTotR,maxTotR]);caxis([minTotR,maxTotR])
    xlim([0.5 1]);
    ylim([0.5 1]);
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    colorbar;
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    subplot(3,3,4)
    surf(ww1,ww2,HSQTFTotRe,'edgecolor','none')
    title(['HYDROSTAR: Re Tot'])
    zlim([minTotR,maxTotR]);caxis([minTotR,maxTotR])
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    colorbar;
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    subplot(3,3,7)
    surf(ww1,ww2,NormalizedAbsDifTOTR,'edgecolor','none')
    title(['Normalized Abs. diff.: Re Tot'])
    %  zlim([0,1]);caxis([0,1])
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    colorbar;
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)

    subplot(3,3,2)
    surf(ww1,ww2,QdatTotI,'edgecolor','none')
    title(['Im Tot'])
    zlim([minTotI,maxTotI]);caxis([minTotI,maxTotI])
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    colorbar;
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    subplot(3,3,5)
    surf(ww1,ww2,HSQTFTotIm,'edgecolor','none')
    title(['- Im Tot'])
    zlim([minTotI,maxTotI]);caxis([minTotI,maxTotI])
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    colorbar;
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    subplot(3,3,8)
    surf(ww1,ww2,NormalizedAbsDifTOTI,'edgecolor','none')
    title(['Im Tot'])
    %  zlim([0,1]);caxis([0,1])
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    colorbar;
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)

    subplot(3,3,3)
    surf(ww1,ww2,QdatTotMod,'edgecolor','none')
    title(['Mod Tot'])
    zlim([minTotMod,maxTotMod]);caxis([minTotMod,maxTotMod])
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    colorbar;
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    subplot(3,3,6)
    surf(ww1,ww2,HSQTFTotMod,'edgecolor','none')
    title(['Mod Tot'])
    zlim([minTotMod,maxTotMod]);caxis([minTotMod,maxTotMod])
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    colorbar;
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    subplot(3,3,9)
    surf(ww1,ww2,NormalizedAbsDifTOTMOD,'edgecolor','none')
    title(['Mod Tot'])
    % zlim([0,1]);caxis([0,1])
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    colorbar;
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)

else
    %% [c,h]=contour plot QTF total
    minQMod=min(min(min(QdatTotMod)));
    maxQMod=max(max(max(QdatTotMod)));

    figure

    if IdContval==1
        [c,h]=contour(ww1,ww2,QdatTotMod,Contval,'ShowText','on','linewidth',1.5);
    else
        [c,h]=contour(ww1,ww2,QdatTotMod,'ShowText','on','linewidth',1.5);
    end
    clabel(c,h,'FontSize',14,'FontWeight','bold')
    title(['NEMOH'])
    xlabel(w1lab);ylabel(w2lab);
    colormap('jet');
    caxis([minQMod, maxQMod])
    %  cb=colorbar;
    %  ylabel(cb,qtflabel,'/\rhog')
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    plot_properties;

    %  cb=colorbar;
    %  ylabel(cb,qtflabel,'/\rhog')
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    plot_properties;

    %% surf plot QTF total
    if  IdMinMaxQTFReIm==0
    minTotI=min(min(QdatTotI));maxTotI=max(max(QdatTotI));
    %     minTotI=-6*1e6;%max(minTotI,min(min(QdatTotR)));
    %     maxTotI=6*1e6;%min(maxTotI,max(max(QdatTotR)));
    else
    minTotI=MinMaxReIm(1);maxTotI=MinMaxReIm(2);
    end
    minTotR=minTotI;%min(min(QdatTotR));
    maxTotR=maxTotI;%max(min(QdatTotR));
    minTotMod=minTotI;%min(min(QdatTotMod));
    maxTotMod=maxTotR;%max(max(QdatTotMod));

    figure;
    subplot(1,3,1)
    surf(ww1,ww2,QdatTotR,'edgecolor','none')
    title(['NEMOH : Re Tot'])
    zlim([minTotR,maxTotR]);caxis([minTotR,maxTotR])
    xlim([0.5 1]);
    ylim([0.5 1]);
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    colorbar;
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)


    subplot(1,3,2)
    surf(ww1,ww2,QdatTotI,'edgecolor','none')
    title(['Im Tot'])
    zlim([minTotI,maxTotI]);caxis([minTotI,maxTotI])
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    colorbar;
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)

    subplot(1,3,3)
    surf(ww1,ww2,QdatTotMod,'edgecolor','none')
    title(['Mod Tot'])
    zlim([minTotMod,maxTotMod]);caxis([minTotMod,maxTotMod])
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    colorbar;
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)

end

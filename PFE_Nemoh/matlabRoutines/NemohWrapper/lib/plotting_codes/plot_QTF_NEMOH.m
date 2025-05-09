    [QdatTotR_surge,QdatTotI_surge,QdatTotMod_surge,w1]=fun_readNEMOHQTFOutput_N...
        ([projdir,filesep,'results',filesep,'QTF',filesep,'OUT_QTF',qtftype,'_N.dat'],...
        NWdatNEM,DOFdatNem,DOFsurge,betaID,NbetaData,qtftype,SwitchBiDir); %data output is already normalized by rho g
    [QdatTotR_heave,QdatTotI_heave,QdatTotMod_heave,w1]=fun_readNEMOHQTFOutput_N...
        ([projdir,filesep,'results',filesep,'QTF',filesep,'OUT_QTF',qtftype,'_N.dat'],...
        NWdatNEM,DOFdatNem,DOFheave,betaID,NbetaData,qtftype,SwitchBiDir); %data output is already normalized by rho g
    [QdatTotR_pitch,QdatTotI_pitch,QdatTotMod_pitch,w1]=fun_readNEMOHQTFOutput_N...
        ([projdir,filesep,'results',filesep,'QTF',filesep,'OUT_QTF',qtftype,'_N.dat'],...
        NWdatNEM,DOFdatNem,DOFpitch,betaID,NbetaData,qtftype,SwitchBiDir); %data output is already normalized by rho g
    
    QdiagR1_surge=diag(QdatTotR_surge,shiftdw1);
    QdiagI1_surge=diag(QdatTotI_surge,shiftdw1);
    QdiagR1_heave=diag(QdatTotR_heave,shiftdw1);
    QdiagI1_heave=diag(QdatTotI_heave,shiftdw1);
    QdiagR1_pitch=diag(QdatTotR_pitch,shiftdw1);
    QdiagI1_pitch=diag(QdatTotI_pitch,shiftdw1);
    
    [ww1,ww2]=meshgrid(w1,w1);
    
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
 
    if strcmpi(qtftype,'P')
        qtflabel='QTF+';
    else
        qtflabel='QTF-';
    end
    freqnorm1=w1(abs(shiftdw1)+1:end);
    dww=w1(2)-w1(1);
    figure;
    subplot(3,1,1)
    plot(freqnorm1,QdiagR1_surge,'-r',freqnorm1,QdiagI1_surge,'--b');
    if ShowLegend==1
        legend('Real','Imag.')
    end
    xlim([min(freqnorm1) max(freqnorm1)])
    if Idwvar==0
        title(['Diag. Surge ', qtflabel,' for \Deltaw=',num2str(abs(shiftdw1)*dww,3)])
    elseif Idwvar==1
        title(['Diag. Surge ', qtflabel,' for \Deltaf=',num2str(abs(shiftdw1)*dww,3)])
    else
        title(['Diag. Surge ', qtflabel,' for \DeltaT=',num2str(abs(shiftdw1)*dww,3)])
    end
    ylabel([qtflabel,'/\rhog'])
    plot_properties;
    subplot(3,1,2)
    plot(freqnorm1,QdiagR1_heave,'-r',freqnorm1,QdiagI1_heave,'--b');
     xlim([min(freqnorm1) max(freqnorm1)])
    ylabel([qtflabel,'/\rhog'])
    if Idwvar==0
        title(['Diag. Heave ', qtflabel,' for \Deltaw=',num2str(abs(shiftdw1)*dww,3)])
    elseif Idwvar==1
        title(['Diag. Heave ', qtflabel,' for \Deltaf=',num2str(abs(shiftdw1)*dww,3)])
    else
        title(['Diag. Heave ', qtflabel,' for \DeltaT=',num2str(abs(shiftdw1)*dww,3)])
    end
    plot_properties;
    subplot(3,1,3)
    plot(freqnorm1,QdiagR1_pitch,'-r',freqnorm1,QdiagI1_pitch,'--b');

    if Idwvar==0
        xlabel('$\omega [rad/s] $','interpreter','latex');
    elseif Idwvar==1
        xlabel('$f[Hz] $','interpreter','latex');
    else
        xlabel('$T [s] $','interpreter','latex');
    end
    xlim([min(freqnorm1) max(freqnorm1)])
    ylabel([qtflabel,'/\rhog'])
    if Idwvar==0
        title(['Diag. Pitch ', qtflabel,' for \Deltaw=',num2str(abs(shiftdw1)*dww,3)])
    elseif Idwvar==1
        title(['Diag. Pitch ', qtflabel,' for \Deltaf=',num2str(abs(shiftdw1)*dww,3)])
    else
        title(['Diag. Pitch ', qtflabel,' for \DeltaT=',num2str(abs(shiftdw1)*dww,3)])
    end
    plot_properties;

    figure;
    subplot(3,3,1)
    surf(ww1,ww2,QdatTotR_surge,'edgecolor','none')
    title(['Surge ',qtflabel,': Real'])
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    cb=colorbar;
    ylabel(cb,[qtflabel,'/\rhog'])
    if Id_fixed_CB_Axis==1
    caxis(minmaxQR_surge)
    end
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    plot_properties;
    subplot(3,3,2)
    surf(ww1,ww2,QdatTotI_surge,'edgecolor','none')
    title(['Imag.'])
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    cb=colorbar;
    ylabel(cb,[qtflabel,'/\rhog'])
    if Id_fixed_CB_Axis==1
    caxis(minmaxQI_surge)
    end
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    plot_properties;
    subplot(3,3,3)
    surf(ww1,ww2,QdatTotMod_surge,'edgecolor','none')
    title(['Mod.'])
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    cb=colorbar;
    ylabel(cb,['|',qtflabel,'|/\rhog'])
    if Id_fixed_CB_Axis==1
    caxis(minmaxQMod_surge)
    end
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    plot_properties;
    subplot(3,3,4)
    surf(ww1,ww2,QdatTotR_heave,'edgecolor','none')
    title(['Heave ',qtflabel,': Real'])
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    cb=colorbar;
    ylabel(cb,[qtflabel,'/\rhog'])
    if Id_fixed_CB_Axis==1
    caxis(minmaxQR_heave)
    end
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    plot_properties;
    subplot(3,3,5)
    surf(ww1,ww2,QdatTotI_heave,'edgecolor','none')
    title(['Imag.'])
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    cb=colorbar;
    ylabel(cb,[qtflabel,'/\rhog'])
    if Id_fixed_CB_Axis==1
    caxis(minmaxQI_heave)
    end
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    plot_properties;
    subplot(3,3,6)
    surf(ww1,ww2,QdatTotMod_heave,'edgecolor','none')
    title(['Mod.'])
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    cb=colorbar;
    ylabel(cb,['|',qtflabel,'|/\rhog'])
    if Id_fixed_CB_Axis==1
    caxis(minmaxQMod_heave)
    end
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    plot_properties;
     subplot(3,3,7)
    surf(ww1,ww2,QdatTotR_pitch,'edgecolor','none')
    title(['Pitch ',qtflabel,': Real'])
    cb=colorbar;
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    ylabel(cb,[qtflabel,'/\rhog'])
    if Id_fixed_CB_Axis==1
    caxis(minmaxQR_pitch)
    end
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    plot_properties;
    subplot(3,3,8)
    surf(ww1,ww2,QdatTotI_pitch,'edgecolor','none')
    title(['Imag.'])
    cb=colorbar;
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    ylabel(cb,[qtflabel,'/\rhog'])
    if Id_fixed_CB_Axis==1
    caxis(minmaxQI_pitch)
    end
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    plot_properties;
    subplot(3,3,9)
    surf(ww1,ww2,QdatTotMod_pitch,'edgecolor','none')
    title(['Mod.'])
    cb=colorbar;
    colormap('jet');xlabel(w1lab);ylabel(w2lab);
    ylabel(cb,['|',qtflabel,'|/\rhog'])
    if Id_fixed_CB_Axis==1
    caxis(minmaxQMod_pitch)
    end
    if Idaxisij==1,axis ij;end
    axis equal
    view(2)
    plot_properties;
hydrostarDir='hydrostar_results_dw0025';
NHeader=22;
ExcF_surge=fun_readHydrostarOutput(['./',hydrostarDir,'/fxf1st.rao'],NHeader);
ExcF_sway=fun_readHydrostarOutput(['./',hydrostarDir,'/fyf1st.rao'],NHeader);
ExcF_heave=fun_readHydrostarOutput(['./',hydrostarDir,'/fzf1st.rao'],NHeader);
ExcF_pitch=fun_readHydrostarOutput(['./',hydrostarDir,'/mxf1st.rao'],NHeader);
ExcF_roll=fun_readHydrostarOutput(['./',hydrostarDir,'/myf1st.rao'],NHeader);
ExcF_yaw=fun_readHydrostarOutput(['./',hydrostarDir,'/mzf1st.rao'],NHeader);
Nw=size(ExcF_surge,1);
Fe_dat=zeros(Nw,13);
Fe_dat(:,[1:2,8])=ExcF_surge(:,1:3);
Fe_dat(:,[3,9])  =ExcF_sway(:,2:3);
Fe_dat(:,[4,10]) =ExcF_heave(:,2:3);
Fe_dat(:,[5,11]) =ExcF_roll(:,2:3);
Fe_dat(:,[6,12]) =ExcF_pitch(:,2:3);
Fe_dat(:,[7,13]) =ExcF_yaw(:,2:3);

Fe_dat(:,8:13)=-Fe_dat(:,8:13);
figure;
subplot(2,1,1)
plot(Fe_dat(:,1),Fe_dat(:,2))
subplot(2,1,2)
plot(Fe_dat(:,1),deg2rad(Fe_dat(:,8)))

save('Fe.dat','Fe_dat','-ascii')
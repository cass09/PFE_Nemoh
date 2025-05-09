hydrostarDir='hydrostar_results_dir30';
NHeader=22;
RAOH_surge=fun_readHydrostarOutput(['./',hydrostarDir,'/surge.rao'],NHeader);
RAOH_sway=fun_readHydrostarOutput(['./',hydrostarDir,'/sway.rao'],NHeader);
RAOH_heave=fun_readHydrostarOutput(['./',hydrostarDir,'/heave.rao'],NHeader);
RAOH_pitch=fun_readHydrostarOutput(['./',hydrostarDir,'/pitch.rao'],NHeader);
RAOH_roll=fun_readHydrostarOutput(['./',hydrostarDir,'/roll.rao'],NHeader);
RAOH_yaw=fun_readHydrostarOutput(['./',hydrostarDir,'/yaw.rao'],NHeader);
Nw=size(RAOH_surge,1);
RAO_dat=zeros(Nw,13);
RAO_dat(:,[1:2,8])=RAOH_surge(:,1:3);
RAO_dat(:,[3,9])  =RAOH_sway(:,2:3);
RAO_dat(:,[4,10]) =RAOH_heave(:,2:3);
RAO_dat(:,[5,11]) =RAOH_roll(:,2:3);
RAO_dat(:,[6,12]) =RAOH_pitch(:,2:3);
RAO_dat(:,[7,13]) =RAOH_yaw(:,2:3);

RAO_dat(:,8:13)=-RAO_dat(:,8:13)+360;
figure;
subplot(2,1,1)
plot(RAO_dat(:,1),RAO_dat(:,3))
subplot(2,1,2)
plot(RAO_dat(:,1),RAO_dat(:,9))

save('RAO.dat','RAO_dat','-ascii')
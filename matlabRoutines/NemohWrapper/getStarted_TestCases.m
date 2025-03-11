%--------------------------------------------------------------------------------------
%
%    Copyright (C) 2022 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   Contributors list:
%   - R. Kurnia
%
%--------------------------------------------------------------------------------------
% This m-file for running the provided NEMOH test-cases
% The input files are already prepared in the TestCases folder
% After running a test case, the results can be verified with the provided reference data
% if Aquaplus data and a tecplot layout .lay are provided in the reference
% folder, copy those files into /result/ then clicking .lay produces
% comparison plots
% if Hydrostar data is provided in the reference folder, uses matlab files in
% postproc_testcase folder for the comparison plots.
%
% User has to specify the TestCases folder in TestCasesDir input below
%-------------------------------------------------------------------------
clc
clear all
close all
[pathstr,~,~] = fileparts(mfilename('fullpath'));
cd (pathstr);
addpath(genpath(pathstr))
%
testcase=1;
TestCasesDir=['.' filesep '..' filesep '..' filesep 'TestCases']; % Needs to be adjusted if you moved the present file
ID_PLOT_RESULTS = 1;
ID_QTF = 0;
ID_HydrosCal = 0;     % A switch,1 computes hydrostatics, inertia, kH.

% Check that Nemoh is available
assert(FindingNemoh(ID_QTF, true))

switch testcase
    case 1     % 1_Cylinder
        projdir=[TestCasesDir,filesep,'1_Cylinder'];
        Nprojdir=1;     %Number of project directory
    case 2     % 2_2Bodies
        projdir=[TestCasesDir,filesep,'2_2Bodies'];
        Nprojdir=1;     %Number of project directory
    case 3     % 3_NonSymmetrical
        projdir=[TestCasesDir,filesep,'3_NonSymmetrical'];
        Nprojdir=1;     %Number of project directory
    case 4     % 4_Postprocessing
        projdir=[TestCasesDir,filesep,'4_Postprocessing'];
        Nprojdir=1;     %Number of project directory
    case 5     % 5_QuickTests
        projdir{1}=[TestCasesDir,filesep,'5_QuickTests',filesep,'1_Sphere'];
        projdir{2}=[TestCasesDir,filesep,'5_QuickTests',filesep,'2_SymmetricSphere'];
        projdir{3}=[TestCasesDir,filesep,'5_QuickTests',filesep,'3_FiniteDepthSphere'];
        projdir{4}=[TestCasesDir,filesep,'5_QuickTests',filesep,'4_SymmetricFiniteDepthSphere'];
        projdir{5}=[TestCasesDir,filesep,'5_QuickTests',filesep,'5_AlienSphere'];
        Nprojdir=5;
        ID_PLOT_RESULTS=0;
    case 6     % 6_box_coarsemesh
        projdir=[TestCasesDir,filesep,'6_box_coarsemesh'];
        Nprojdir=1;     %Number of project directory
    case 7     % 7_Solvers_Check_OC3
        projdir{1}=[TestCasesDir,filesep,'7_Solvers_Check_OC3',filesep,'mesh978_GaussElim'];
        projdir{2}=[TestCasesDir,filesep,'7_Solvers_Check_OC3',filesep,'mesh978_GMRES'];
        projdir{3}=[TestCasesDir,filesep,'7_Solvers_Check_OC3',filesep,'mesh978_LUdecomp'];
        projdir{4}=[TestCasesDir,filesep,'7_Solvers_Check_OC3',filesep,'mesh4842_GaussElim'];
        projdir{5}=[TestCasesDir,filesep,'7_Solvers_Check_OC3',filesep,'mesh4842_GMRES'];
        projdir{6}=[TestCasesDir,filesep,'7_Solvers_Check_OC3',filesep,'mesh4842_LUdecomp'];
        Nprojdir=6;     %Number of project directory
        ID_PLOT_RESULTS=0;
    case 8      % 8a_Cylinder_irregfreq
        projdir{1}=[TestCasesDir,filesep,'8a_Cylinder_irregfreq',filesep,'mesh634'];
        projdir{2}=[TestCasesDir,filesep,'8a_Cylinder_irregfreq',filesep,'mesh634_nolid'];
        Nprojdir=2;     %Number of project directory
    case 9      % 8b_QTF_Cylinder
        projdir=[TestCasesDir,filesep,'8b_QTF_Cylinder',filesep,'mesh634_floating'];
        Nprojdir=1;     %Number of project directory
        ID_HydrosCal=1; % if error due to mesh.cal please copy the Mesh.cal
        % in the same folder as this m-file.
        % The error due to the executables are not updated
        % yet.
        ID_QTF=1;
    case 10      % 9_QTF_OC4_Semisubmersible
        projdir=[TestCasesDir,filesep,'9_QTF_OC4_Semisubmersible',filesep,'mesh2200_Floating'];
        Nprojdir=1;     %Number of project directory
        ID_QTF=1;
    case 11      % 10a_QTF_SOFTWIND
        projdir=[TestCasesDir,filesep,'10_QTF_SOFTWIND',filesep,'mesh_1872half_QTF_Floating'];
        Nprojdir=1;     %Number of project directory
        ID_QTF=1;
    case 12      % 10b_QTF_SOFTWIND
        projdir=[TestCasesDir,filesep,'10_QTF_SOFTWIND',filesep,'mesh_1872half_nolid_QTF_Floating_FS'];
        Nprojdir=1;     %Number of project directory
        ID_QTF=1;
    case 13      % 11_QTF_OC3_Hywind
        projdir=[TestCasesDir,filesep,'11_QTF_OC3_Hywind',filesep,'mesh4842_Floating'];
        Nprojdir=1;     %Number of project directory
        ID_HydrosCal=1; % if error due to mesh.cal please copy the Mesh.cal
        % in the same folder as this m-file.
        % The error due to the executables are not updated
        % yet.
        ID_QTF=1;
end


%-------Launch Calculation------------
if Nprojdir==1
    [Idw,w,A,B,Fe]=Nemoh(projdir,ID_HydrosCal,ID_QTF); % Call the function Nemoh.m
else
    for iproj=1:Nprojdir
        [Idw,w,A,B,Fe]=Nemoh(projdir{iproj},ID_HydrosCal,ID_QTF); % Call the function Nemoh.m
    end
end
%% --- Computes QTFs --------------------
if ID_QTF==1
    if Nprojdir==1
        NemohQTF(projdir); % Call the function NemohQTF.m
    else
        for iproj=1:Nprojdir
            NemohQTF(projdir{iproj}); % Call the function NemohQTF.m
        end
    end
end
%% ---------------------------------------
if ID_PLOT_RESULTS==1 %Plot NEMOH1
    if Idw==1
        xlab='\omega [rad/s]';
    elseif Idw==2
        xlab='f [Hz]';
    elseif Idw==3
        xlab='T [s]';
    end
    
    figure
    subplot(3,2,1)
    a(1,:)=A(1,1,:);
    b(1,:)=B(1,1,:);
    plot(w,a,'b-+',w,b,'r-+')
    grid ON
    ylabel('Coeffs' )
    legend('Added mass (A)','Damping (B)')
    title('Surge')
    plot_properties;
    
    subplot(3,2,2)
    plot(w,abs(Fe(:,1)),'g-+')
    grid ON
    ylabel('F (N)' )
    title('Surge')
    plot_properties;
    
    subplot(3,2,3)
    a(1,:)=A(3,3,:);
    b(1,:)=B(3,3,:);
    plot(w,a,'b-+',w,b,'r-+')
    grid ON
    ylabel('Coeffs' )
    legend('Added mass (A)','Damping (B)')
    title('Heave')
    plot_properties;
    
    subplot(3,2,4)
    plot(w,abs(Fe(:,3)),'g-+')
    grid ON
    ylabel('F(N)' )
    title('Heave')
    plot_properties;
    
    subplot(3,2,5)
    a(1,:)=A(5,5,:);
    b(1,:)=B(5,5,:);
    plot(w,a,'b-+',w,b,'r-+')
    xlabel(xlab)
    ylabel('Coeffs' )
    grid ON
    legend('Added mass (A)','Damping (B)')
    title('Pitch')
    plot_properties;
    
    subplot(3,2,6)
    plot(w,abs(Fe(:,5)),'g-+')
    xlabel(xlab)
    grid ON
    ylabel('F(N)' )
    title('Pitch')
    plot_properties;
    
    if ID_QTF==1
        qtftype='M';
        Idwvar=0;%0=w rad/s,1=f 1/s,2=T s; %freq type
        SwitchBiDir=0;
        betaID=[0 0]; %[beta1 beta2]
        shiftdw1=-1;
        ShowLegend=1;
        Idaxisij=0;
        NWdatNEM=65;%specify number of freq as in the QTF output
        DOFdatNem=[1 1 1 1 1 1];% the QTF data availability in each mode, 1 means data exist
        % if for index 1 to 6, all the values are 1, it
        % means data available for all DOF
        DOFsurge=1;DOFheave=3;DOFpitch=5;
        NbetaData=1;
        
        Id_fixed_CB_Axis=0; %switch to fixed colorbar limit, 
                            % if 1 the minmaxQ values have to be specified
        minmaxQR_surge=[-10 30]; minmaxQI_surge=[-60 60];minmaxQMod_surge=[0 60];
        minmaxQR_heave=[-5 5];   minmaxQI_heave=[-5 5];  minmaxQMod_heave=[0 5];
        minmaxQR_pitch=[-50 150]; minmaxQI_pitch=[-50 150];minmaxQMod_pitch=[0 150];

        plot_QTF_NEMOH;
    end
end
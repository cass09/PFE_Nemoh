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
%--------------------------------------------------------------------------------------
%
% --> function []=NemohQTF(projdir)
%
% Purpose: NEMOH Matlab wrapper for computing QTFs
%
% Inputs :
% - projdir: Project directory
%
%
function NemohQTF(projdir)

system(['mkdir ',projdir,filesep,'QTFPreprocOut']);
system(['mkdir ',projdir,filesep,'results',filesep,'QTF']);

% Calcul des coefficients hydrodynamiques
fprintf('\n------ Starting QTF preproc ----------- \n');
system(['QTFpreProc ',projdir]);
fprintf('------ Computing QTFs ------------- \n');
system(['QTFsolver ',projdir]);
fprintf('------ Postprocessing results --- \n');
system(['QTFpostProc ',projdir]);

end

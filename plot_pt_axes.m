% PLOT_PT_AXES Plot PT-axes into polarplot diagram.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Plot PT-axes into polarplot diagram
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Miroslav HALLO
% Charles University in Prague, Faculty of Mathematics and Physics
% E-mail: hallo@karel.troja.mff.cuni.cz
% Revision 2018/07: The first version of the function
% Revision 2018/12: Enhanced version
% Revision 2026/04: New Matlab version
% Tested in Matlab R2025b
% Method:
% Hallo,M., Oprsal,I., Asano,K., Gallovic,F. (2019): Seismotectonics of the 2018
%      Northern Osaka M6.1 earthquake and its aftershocks: joint movements
%      on strike-slip and reverse faults in inland Japan, Earth,
%      Planets and Space, 71:34. https://doi.org/10.1186/s40623-019-1016-8
%
% Copyright (C) 2018 Miroslav Hallo
%
% This program is published under the GNU General Public License (GNU GPL).
%
% This program is free software: you can modify it and/or redistribute it
% or any derivative version under the terms of the GNU General Public
% License as published by the Free Software Foundation, either version 3
% of the License, or (at your option) any later version.
%
% This code is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
% and don't remove their names from the code.
%
% You should have received copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INIT:
close all;
clearvars;
projRoot = fileparts(which(mfilename));
addpath(fullfile(projRoot, 'lib'));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:

% Input text file with Strike, Dip, Rake of double-couple focal mechanisms [deg]
sdr_file = 'example_sdr.txt';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create directory for results
resDir = fullfile(projRoot, 'results');
if ~exist(resDir, 'dir')
    mkdir(resDir);
end

% Prepare timestamp
timestamp = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
outfile = [char(timestamp),'_pt_axes'];

% Read file with Strike, Dip, Rake
try
    fid = fopen(sdr_file,'r');
    sdrlist = textscan(fid,'%f %f %f %[^\n]', 'CommentStyle', '#');
    fclose(fid);
catch
    sdrlist = [];
end
total = length(sdrlist{1,1});
sdr = [sdrlist{1,1} sdrlist{1,2} sdrlist{1,3}];

% Prepare Theta and Rho of all PT-axes
[P_polarThe,P_polarRho,T_polarThe,T_polarRho] = mech_pt(sdr(:,1),sdr(:,2),sdr(:,3));

% Prepare polar axes
fang = 0 : 0.1 : 360;
frad = fang .* pi/180;

% Plot figure
figure('color','white');

polarplot(frad, ones(1,length(fang)), 'k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
hold on;
polarplot(P_polarThe, P_polarRho, 'ro', 'MarkerSize', 5, 'LineWidth', 1.0);
polarplot(T_polarThe, T_polarRho,'b+', 'MarkerSize', 5, 'LineWidth', 1.0);

ax = gca;
ax.ThetaTick = 0:22.5:360;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.ThetaTickLabels = {'N','','NE','','E','','SE','','S','','SW','','W','','NW',''};
ax.RTickLabel = [];
ax.GridColor = [0.9 0.9 0.9];
ax.GridAlpha = 1;
hold off

legend('P-axis','T-axis','Location','northeastoutside')

% Save image
exportgraphics(gcf, fullfile(resDir, [outfile,'.png']), 'Resolution', 300);
exportgraphics(gcf, fullfile(resDir, [outfile,'.pdf']), 'ContentType', 'vector');

fprintf('Figures successfully saved as: %s.png and %s.pdf\n', outfile, outfile);


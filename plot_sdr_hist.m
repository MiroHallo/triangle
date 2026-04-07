% PLOT_SDR_HIST Plot polar histograms of strike and dip angles.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot polar histograms of strike and dip angles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Miroslav HALLO, Ivo OPRSAL
% Charles University in Prague, Faculty of Mathematics and Physics
% E-mail: hallo@karel.troja.mff.cuni.cz
% Revision 2018/07: The first version of the function
% Revision 2018/08: Enhanced version
% Revision 2026/04: New Matlab version
% Tested in Matlab R2025b
% Method:
% Hallo,M., Oprsal,I., Asano,K., Gallovic,F. (2019): Seismotectonics of the 2018
%      Northern Osaka M6.1 earthquake and its aftershocks: joint movements
%      on strike-slip and reverse faults in inland Japan, Earth,
%      Planets and Space, 71:34. https://doi.org/10.1186/s40623-019-1016-8
%
% Copyright (C) 2018 Miroslav Hallo and Ivo Oprsal
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT:

% Input text file with Strike, Dip, Rake of double-couple focal mechanisms [deg]
sdr_file = 'example_sdr.txt';

% Dip-slice style
% 0 = 0-90,   1 = W-E,   2 = N-S,   3 = NW-SE,   4 = SW-NE
% Note: Select 0 or the perpendicular direction to the main strike azimuth
dipSlice = 3;

% Binsize for the polar histogram (deg)
% Number of bins is rounded to neares integer nbins = round(360/anglestep)
anglestep = 15;

% Threshold to filter out events with small dip angle (deg)
dip_th = 5;

% Set RGB color
coloF = [0.30, 0.30, 0.30];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create directory for results
resDir = fullfile(projRoot, 'results');
if ~exist(resDir, 'dir')
    mkdir(resDir);
end

% Prepare timestamp
timestamp = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
outfile = [char(timestamp),'_sd_hist'];

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

% Prepare polar axis
deg2rad = pi/180;
nbins = round(360/anglestep);
anglestep = 360/nbins;
bins = 0:anglestep:360;
polar_data = [sdr(:,2),sdr(:,1),sdr(:,3)];
polar_data = polar_data(polar_data(:,1)>dip_th,:);

if dipSlice==1
    TickL = 'E'; TickR = 'W';
    for i=1:length(polar_data(:,1))
        if (polar_data(i,2)>90) && (polar_data(i,2)<270)
            polar_data(i,1) = 180 - polar_data(i,1);
        end
    end
elseif dipSlice==2
    TickL = 'S'; TickR = 'N';
    for i=1:length(polar_data(:,1))
        if (polar_data(i,2)>180)
            polar_data(i,1) = 180 - polar_data(i,1);
        end
    end
elseif dipSlice==3
    TickL = 'SE'; TickR = 'NW';
    for i=1:length(polar_data(:,1))
        if (polar_data(i,2)>135) && (polar_data(i,2)<315)
            polar_data(i,1) = 180 - polar_data(i,1);
        end
    end
elseif dipSlice==4
    TickL = 'NE'; TickR = 'SW';
    for i=1:length(polar_data(:,1))
        if (polar_data(i,2)>45) && (polar_data(i,2)<225)
            polar_data(i,1) = 180 - polar_data(i,1);
        end
    end
end

% Flipping and mirroring strikes
polar_data(polar_data(:,2)>=180,2) = polar_data(polar_data(:,2)>=180,2) - 180;
tmp = polar_data;
tmp(:,2) = tmp(:,2) + 180;
polar_data = [polar_data' tmp']';
clearvars tmp;

% Plot azimuths
figure('color','white');
subplot(1,2,1)
h = polarhistogram(polar_data(:,2)*deg2rad, bins*deg2rad, 'FaceColor', coloF, ...
    'EdgeColor','k','FaceAlpha',1,'DisplayStyle','bar','Normalization','probability');
maxNum = max(h.Values);
rTicks = [maxNum/4 maxNum/2 3*maxNum/4 maxNum];
ax = gca;
ax.RLim = [0 maxNum*1.1];
ax.RTick = rTicks;
ax.ThetaTick=[0 15 30 45 60 75 90 105 120 135 150 165 180 195 210 225 240 255 270 285 300 315 330 345];
ax.ThetaTickLabels = {'N' '' '30' '' '60' '' 'E' '' '120' '' '150' '' 'S' '' '210' '' '240' '' 'W' '' '300' '' '330' ''};
ax.GridAlpha = 1;
ax.GridColor = [0.9 0.9 0.9];
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
ax.RTick = rTicks;
ax.RTickLabel = [];
title(['Strike azimuths (',num2str(length(polar_data(:,1))/2),' nodal planes)'],...
    'FontWeight', 'normal')

% Plot dips
subplot(1,2,2)
h = polarhistogram(polar_data(:,1)*deg2rad, bins*deg2rad, 'FaceColor', coloF, ...
    'EdgeColor','k','FaceAlpha',1,'DisplayStyle','bar','Normalization','probability');
maxNum = max(h.Values);
rTicks = [maxNum/4 maxNum/2 3*maxNum/4 maxNum];
ay = gca;
ay.RLim = [0 maxNum*1.1];
if dipSlice~=0
    ay.ThetaLim=[0 180];
    ay.ThetaTick=[0 15 30 45 60 75 90 105 120 135 150 165 180];
    ay.ThetaTickLabels = {TickL,'','30','','60','','DOWN','','60','','30','',TickR};
else
    ay.ThetaLim=[0 90];
    ay.ThetaTick=[0 15 30 45 60 75 90];
    ay.ThetaTickLabels = {'0','','30','','60','','DOWN'};
end
ay.GridAlpha = 1;
ay.GridColor = [0.9 0.9 0.9];
ay.ThetaDir = 'clockwise';
ay.ThetaZeroLocation = 'right';
ay.RTick = rTicks;
ay.RTickLabel = [];
title(['Dip angles (',num2str(length(polar_data(:,1))/2),' nodal planes)'],...
    'FontWeight', 'normal')

% Save image
exportgraphics(gcf, fullfile(resDir, [outfile,'.png']), 'Resolution', 300);
exportgraphics(gcf, fullfile(resDir, [outfile,'.pdf']), 'ContentType', 'vector');

fprintf('Figures successfully saved as: %s.png and %s.pdf\n', outfile, outfile);


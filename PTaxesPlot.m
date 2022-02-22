%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Plot PT-axes into polarplot diagram 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Coded for the purpose of paper Hallo et al. (2019)
% Hallo,M., Oprsal,I., Asano,K., Gallovic,F. (2019): Seismotectonics of the 2018
%      Northern Osaka M6.1 earthquake and its aftershocks: joint
%      movements on strike-slip and reverse faults in inland Japan, Earth,
%      Planets and Space, 71:34.
%
% Code author: Miroslav Hallo
% Charles University in Prague, Faculty of Mathematics and Physics
% Web: http://geo.mff.cuni.cz/~hallo/
% E-mail: hallo@karel.troja.mff.cuni.cz
% Revision 7/2018: The first version of the function.
% Revision 12/2018: Enhanced version.
%
% Copyright (C) 2018  Miroslav Hallo
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

close all;
clear all;
addpath([pwd,'/MATLAB']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT
% List of Strike / Dip / Rakes:
SDR = [
254	83	-163
37	59	84
25	52	86
40	45	95
236	72	-158
247	65	-168
26	59	77
6	32	78
241	81	-134
49	46	97
251	50	160
230	82	-137
241	86	-160
351	54	44
352	22	67
247	63	-167
33	46	91
53	65	157
250	84	-163
71	51	132
];

% END INPUT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[P_polarThe,P_polarRho,T_polarThe,T_polarRho] = mechPT(SDR(:,1),SDR(:,2),SDR(:,3));
N = length(SDR(:,1));

figure('color','white');
try % MATLAB version > 2016 (new)
    fang = 0:0.1:360;
    polarplot(fang*pi/180,ones(1,length(fang)),'k','LineWidth',1.0,'HandleVisibility','off');
    hold on;
    polarplot(P_polarThe,P_polarRho,'ro','MarkerSize',5,'LineWidth',1.0);
    polarplot(T_polarThe,T_polarRho,'b+','MarkerSize',5,'LineWidth',1.0);
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.ThetaTick = 0:45:360;
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
    ax.ThetaTickLabels = {'N','NE','E','SE','S','SW','W','NW'};
    ax.FontWeight = 'bold';
    ax.RTickLabel = [];
    hold off;
    legend('P-axis','T-axis','Location','northeastoutside')
catch % MATLAB version < 2016 (old)
    fang = 0:0.1:360;
    hp = mmpolar(fang*pi/180,ones(1,length(fang)),'k',P_polarThe,P_polarRho,'ro',T_polarThe,T_polarRho,'b+', ...
        'Style','Compass', 'TTickDelta',45,'TGridLineStyle','-','TGridColor',[0.8 0.8 0.8], ...
        'TTickLabel',{'N' 'NE' 'E' 'SE' 'S' 'SW' 'W' 'NW'},'RLimit',[0,1], ...
        'RGridLineStyle','-','RGridColor',[0.8 0.8 0.8],'RTickLabelVisible','off', ...
        'FontWeight','bold');
    legend(hp(2:3),'P-axis','T-axis','Location','northeastoutside')
end





return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Classification of source focal mechanisms and plot the triangle diagram 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Method by Frohlich (1992)
% Frohlich,C. (1992): Triangle diagrams: ternary graphs to display similarity
%      and diversity of earthquake focal mechanisms, Physics of the Earth and 
%      Planetary Interiors, 75, 193-198.
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
% Revision 2/2019: Plot the class boundaries and additional "Odd" class legends.
%
% Copyright (C) 2018,2019  Miroslav Hallo
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

% Color RGB codes (odd, strike-slip, normal, reverse)
coloF = [
0.2 0.2 0.2
1 0 0
0 1 0
0 0 1
];

% Plot grid (strike-slip, normal, reverse)
pGrid = 1;

% Plot secondary grid (90, 80, ... deg)
pSecGrid = 1;

% Plot additional legend of odd focal mechanisms
pAddOdd = 1;

% END INPUT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dN = 35.26; 
[mClassAll,dPAll,dTAll,dBAll] = mechClass(SDR(:,1),SDR(:,2),SDR(:,3));

figure('color','w');
hold on

% Plot Border
dP = 90; dT = 0; dB = 0;
z = atand(sind(dT)/sind(dP))-45;
h1 = (cosd(dB)*sind(z)) / (sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
v1 = (cosd(dN)*sind(dB)-sind(dN)*cosd(dB)*cosd(z))/(sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));

dP = 0; dT = 90; dB = 0;
z = atand(sind(dT)/sind(dP))-45;
h2 = (cosd(dB)*sind(z)) / (sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
v2 = (cosd(dN)*sind(dB)-sind(dN)*cosd(dB)*cosd(z))/(sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));

dP = 0; dT = 0; dB = 90;
z = atand(1)-45;
h3 = (cosd(dB)*sind(z)) / (sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
v3 = (cosd(dN)*sind(dB)-sind(dN)*cosd(dB)*cosd(z))/(sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
plot([h1 h2 h3 h1],[v1 v2 v3 v1],'-k')

% Plot Secondary Grid
if pSecGrid
    Np = 100;
    for dd = 20:10:80
        h1_grid = zeros(1,Np);
        v1_grid = zeros(1,Np);
        for i=1:Np
            dP = dd;
            dT = (i-1)*((90-dd)/(Np-1));
            dB = asind(sqrt(1 - (sind(dP)^2 + sind(dT)^2)));
            z = atand(sind(dT)/sind(dP))-45;
            h1_grid(i) = (cosd(dB)*sind(z)) / (sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
            v1_grid(i) = (cosd(dN)*sind(dB)-sind(dN)*cosd(dB)*cosd(z))/(sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
        end
        plot(h1_grid,v1_grid,'-','Color',[0.6 0.6 0.6],'LineWidth',0.5)
        text(h1_grid(end),v1_grid(end)-0.08,num2str(dd),'VerticalAlignment','top','HorizontalAlignment','center','Color',[0.6 0.6 0.6])
        
        h2_grid = zeros(1,Np);
        v2_grid = zeros(1,Np);
        for i=1:Np
            dB = dd;
            dT = (i-1)*((90-dd)/(Np-1));
            dP = asind(sqrt(1 - (sind(dB)^2 + sind(dT)^2)));
            z = atand(sind(dT)/sind(dP))-45;
            h2_grid(i) = (cosd(dB)*sind(z)) / (sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
            v2_grid(i) = (cosd(dN)*sind(dB)-sind(dN)*cosd(dB)*cosd(z))/(sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
        
        end
        [h2_grid,ind] = sort(h2_grid);
        plot(h2_grid,v2_grid(ind),'-','Color',[0.6 0.6 0.6],'LineWidth',0.5)
        text(h2_grid(1)-0.1,v2_grid(1),num2str(dd),'VerticalAlignment','bottom','HorizontalAlignment','right','Color',[0.6 0.6 0.6])
        
        h3_grid = zeros(1,Np);
        v3_grid = zeros(1,Np);
        for i=1:Np
            dT = dd;
            dP = (i-1)*((90-dd)/(Np-1));
            dB = asind(sqrt(1 - (sind(dP)^2 + sind(dT)^2)));
            z = atand(sind(dT)/sind(dP))-45;
            h3_grid(i) = (cosd(dB)*sind(z)) / (sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
            v3_grid(i) = (cosd(dN)*sind(dB)-sind(dN)*cosd(dB)*cosd(z))/(sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
        end
        plot(h3_grid,v3_grid,'-','Color',[0.6 0.6 0.6],'LineWidth',0.5)
        text(h3_grid(1)+0.1,v3_grid(1),num2str(dd),'VerticalAlignment','bottom','HorizontalAlignment','left','Color',[0.6 0.6 0.6])
    end
    
    text(0,-1.08,'\delta_P (deg)','VerticalAlignment','bottom','HorizontalAlignment','center','Color',[0.6 0.6 0.6])
    text(-0.9,0.5,'\delta_B (deg)','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',60,'Color',[0.6 0.6 0.6])
    text(0.9,0.5,'\delta_T (deg)','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',-60,'Color',[0.6 0.6 0.6])
    
end

% Plot Major Grid
if pGrid
    Np = 100;
    plot(0,0,'+','Color',[0.6 0.6 0.6])
    
    h1_grid = zeros(1,Np);
    v1_grid = zeros(1,Np);
    for i=1:Np
        dP = 60;
        dT = (i-1)*(30/(Np-1));
        dB = asind(sqrt(1 - (sind(dP)^2 + sind(dT)^2)));
        z = atand(sind(dT)/sind(dP))-45;
        h1_grid(i) = (cosd(dB)*sind(z)) / (sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
        v1_grid(i) = (cosd(dN)*sind(dB)-sind(dN)*cosd(dB)*cosd(z))/(sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
    end
    plot(h1_grid,v1_grid,'-','Color',[0.6 0.6 0.6],'LineWidth',2)
    
    h2_grid = zeros(1,Np);
    v2_grid = zeros(1,Np);
    for i=1:Np
        dB = 60;
        dT = (i-1)*(30/(Np-1));
        dP = asind(sqrt(1 - (sind(dB)^2 + sind(dT)^2)));
        z = atand(sind(dT)/sind(dP))-45;
        h2_grid(i) = (cosd(dB)*sind(z)) / (sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
        v2_grid(i) = (cosd(dN)*sind(dB)-sind(dN)*cosd(dB)*cosd(z))/(sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
    end
    [h2_grid,ind] = sort(h2_grid);
    plot(h2_grid,v2_grid(ind),'-','Color',[0.6 0.6 0.6],'LineWidth',2)
    
    h3_grid = zeros(1,Np);
    v3_grid = zeros(1,Np);
    for i=1:Np
        dT = 50;
        dP = (i-1)*(40/(Np-1));
        dB = asind(sqrt(1 - (sind(dP)^2 + sind(dT)^2)));
        z = atand(sind(dT)/sind(dP))-45;
        h3_grid(i) = (cosd(dB)*sind(z)) / (sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
        v3_grid(i) = (cosd(dN)*sind(dB)-sind(dN)*cosd(dB)*cosd(z))/(sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
    end
    plot(h3_grid,v3_grid,'-','Color',[0.6 0.6 0.6],'LineWidth',2)
    
end

% Plot Events
for i=1:length(mClassAll)
    dP = dPAll(i);
    dT = dTAll(i);
    dB = dBAll(i);
    z = atand(sind(dT)/sind(dP))-45;
    
    h = (cosd(dB)*sind(z)) / (sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
    v = (cosd(dN)*sind(dB)-sind(dN)*cosd(dB)*cosd(z))/(sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
    plot(h,v,'o','MarkerFaceColor',coloF(mClassAll(i)+1,:),'MarkerEdgeColor',coloF(mClassAll(i)+1,:),'MarkerSize',7)
end

% Plot Focal Mechanisms
bb([0,45,-90], h1, v1, 0.12, 0, coloF(3,:))
bb([0,45,90], h2, v2, 0.12, 0, coloF(4,:))
bb([45,90,180], h3, v3, 0.12, 0, coloF(2,:))

% Plot Additional Focal Mechanisms
if pAddOdd
    SDRodd = [0,90,-90];
    [~,dP,dT,dB] = mechClass(SDRodd(1),SDRodd(2),SDRodd(3));
    z = atand(sind(dT)/sind(dP))-45;
    hodd = (cosd(dB)*sind(z)) / (sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
    vodd = (cosd(dN)*sind(dB)-sind(dN)*cosd(dB)*cosd(z))/(sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
    bb(SDRodd, hodd, vodd, 0.12, 0, coloF(1,:))
    
    SDRodd = [144,60,35];
    [~,dP,dT,dB] = mechClass(SDRodd(1),SDRodd(2),SDRodd(3));
    z = atand(sind(dT)/sind(dP))-45;
    hodd = (cosd(dB)*sind(z)) / (sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
    vodd = (cosd(dN)*sind(dB)-sind(dN)*cosd(dB)*cosd(z))/(sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
    bb(SDRodd, hodd, vodd, 0.12, 0, coloF(1,:))
    
    SDRodd = [216,60,-35];
    [~,dP,dT,dB] = mechClass(SDRodd(1),SDRodd(2),SDRodd(3));
    z = atand(sind(dT)/sind(dP))-45;
    hodd = (cosd(dB)*sind(z)) / (sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
    vodd = (cosd(dN)*sind(dB)-sind(dN)*cosd(dB)*cosd(z))/(sind(dN)*sind(dB)+cosd(dN)*cosd(dB)*cosd(z));
    bb(SDRodd, hodd, vodd, 0.12, 0, coloF(1,:))
end

hold off
axis equal
axis off




return
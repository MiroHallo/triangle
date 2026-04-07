function [P_polarThe,P_polarRho,T_polarThe,T_polarRho] = mech_pt(strike,dip,rake)
% MECH_PT Returns PT-axes in polarplot Theta and Rho.
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
%
% INPUT:
% strike - Strike angle of one nodal plane [deg]
% dip - Dip angle of one nodal plane [deg]
% rake - Rake angle of one nodal plane [deg]
% 
% OUTPUT:
% P_polarThe, P_polarRho - polarplot Theta and Rho of P-axis
% T_polarThe, T_polarRho - polarplot Theta and Rho of T-axis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normals and slip vectors
n0(:,1) = -sind(dip).*sind(strike);
n0(:,2) =  sind(dip).*cosd(strike);
n0(:,3) = -cosd(dip);
u0(:,1) =  cosd(rake).*cosd(strike) + cosd(dip).*sind(rake).*sind(strike);
u0(:,2) =  cosd(rake).*sind(strike) - cosd(dip).*sind(rake).*cosd(strike);
u0(:,3) = -sind(rake).*sind(dip);

%-------------------------------------------------------------------
% PT-axes
P_osa = (n0-u0) ./ repmat(sqrt(sum((n0-u0).^2,2)),1,3);
T_osa = (n0+u0) ./ repmat(sqrt(sum((n0+u0).^2,2)),1,3);
P_osa(P_osa(:,3)>0,:) = -P_osa(P_osa(:,3)>0,:);
T_osa(T_osa(:,3)>0,:) = -T_osa(T_osa(:,3)>0,:);

%-------------------------------------------------------------------
% Compute all angles
N = length(strike);
P_azimuth = zeros(N,1);
T_azimuth = zeros(N,1);
P_theta = zeros(N,1);
T_theta = zeros(N,1);
P_x = zeros(N,1);
P_y = zeros(N,1);
T_x = zeros(N,1);
T_y = zeros(N,1);

for i=1:N
    % Get azimuths and dip angles
    P_azimuth(i) = atan2d(P_osa(i,1),P_osa(i,2));
    P_theta(i) = acosd(abs(P_osa(i,3)));
    
    T_azimuth(i) = atan2d(T_osa(i,1),T_osa(i,2));
    T_theta(i) = acosd(abs(T_osa(i,3)));
    
    % Convert to projection
    P_x(i) = -sqrt(2.) * sind(P_theta(i)/2) .* sind(P_azimuth(i));
    P_y(i) = -sqrt(2.) * sind(P_theta(i)/2) .* cosd(P_azimuth(i));
    
    T_x(i) = -sqrt(2.) * sind(T_theta(i)/2) .* sind(T_azimuth(i));
    T_y(i) = -sqrt(2.) * sind(T_theta(i)/2) .* cosd(T_azimuth(i));
end

% Convert to polar coordinates and plot
[P_polarThe, P_polarRho] = cart2pol(P_x, P_y);
[T_polarThe, T_polarRho] = cart2pol(T_x, T_y);

end

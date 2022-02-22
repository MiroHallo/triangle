%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Returns PT-axes in polarplot Theta and Rho 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Coded for the purpose of paper Hallo et al. (2019)
% Hallo,M., Oprsal,I., Asano,K., Gallovic,F. (2019): Seismotectonics of the 2018
%      Northern Osaka M6.1 earthquake and its aftershocks: joint
%      movements on strike-slip and reverse faults in inland Japan, Earth,
%      Planets and Space, submitted.
%
% Code author: Miroslav Hallo
% Charles University in Prague, Faculty of Mathematics and Physics
% Web: http://geo.mff.cuni.cz/~hallo/
% E-mail: hallo@karel.troja.mff.cuni.cz
% Revision 7/2018: The first version of the function.
% Revision 12/2018: Enhanced version.
%
% This code is published under the GNU General Public License. To any
% licensee is given permission to modify the work, as well as to copy
% and redistribute the work or any derivative version. Still we would
% like to kindly ask you to acknowledge the authors and don't remove
% their names from the code. This code is distributed in the hope
% that it will be useful, but WITHOUT ANY WARRANTY.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P_polarThe,P_polarRho,T_polarThe,T_polarRho] = mechPT(strike,dip,rake)
% Returns:
% P_polarThe,P_polarRho - polarplot Theta and Rho of P-axis
% T_polarThe,T_polarRho - polarplot Theta and Rho of T-axis

%--------------------------------------------------------------------------
% Normals and slip vectors
n0(:,1) = -sind(dip).*sind(strike);
n0(:,2) =  sind(dip).*cosd(strike);
n0(:,3) = -cosd(dip);
u0(:,1) =  cosd(rake).*cosd(strike) + cosd(dip).*sind(rake).*sind(strike);
u0(:,2) =  cosd(rake).*sind(strike) - cosd(dip).*sind(rake).*cosd(strike);
u0(:,3) = -sind(rake).*sind(dip);

%--------------------------------------------------------------------------
% PT-axes
P_osa = (n0-u0)./repmat(sqrt(sum((n0-u0).^2,2)),1,3);
T_osa = (n0+u0)./repmat(sqrt(sum((n0+u0).^2,2)),1,3);
P_osa(P_osa(:,3)>0,:) = -P_osa(P_osa(:,3)>0,:);
T_osa(T_osa(:,3)>0,:) = -T_osa(T_osa(:,3)>0,:);

%--------------------------------------------------------------------------
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
    P_x(i) = -sqrt(2.)*sind(P_theta(i)/2).*sind(P_azimuth(i));
    P_y(i) = -sqrt(2.)*sind(P_theta(i)/2).*cosd(P_azimuth(i));
    
    T_x(i) = -sqrt(2.)*sind(T_theta(i)/2).*sind(T_azimuth(i));
    T_y(i) = -sqrt(2.)*sind(T_theta(i)/2).*cosd(T_azimuth(i));
end

% Convert to polar coordinates and plot
[P_polarThe,P_polarRho] = cart2pol(P_x,P_y);
[T_polarThe,T_polarRho] = cart2pol(T_x,T_y);

end

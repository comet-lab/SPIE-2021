function [h,u] = calc_config(L, R, n, ro, ri, w)
%% Calculates the configuration of a tube given desired characteristics
%  [h,u] = calc_config(L, R, n, ro, ri) returns h and u for a tube 
%  [h,u] = calc_config(L, R, n) returns h and u for a tube, assuming a
%  default tube
%   L [mm] = steerable section length (max arc length)
%   R [mm] = max radius of curvature
%   n = number of notches
%   *optional* 
%   ro [mm] = outer radius of tube
%   ri [mm] = inner radius of tube
%   w [mm] = depth of notches
%
%   Last Revision: 5/27/2020

%% define optional variables

if ~exist('ro', 'var')
    % set default outer radius
    ro = 1.6/2;
end
if ~exist('ri', 'var')
    % set default inner radius
    ri = 1.4/2;
end
if ~exist('w', 'var')
    % set default notch depth percentage
    w = 1.4; 
end
%% calculate ybar
phi_o = 2*acos((w-ro)/ro);
phi_i = 2*acos((w-ro)/ri);

Ao = (ro^2*(phi_o - sin(phi_o)))/2;
Ai = (ri^2*(phi_i - sin(phi_i)))/2;

ybar_o = (4*ro*sin(phi_o/2)^3)/(3*(phi_o - sin(phi_o)));
ybar_i = (4*ri*sin(phi_i/2)^3)/(3*(phi_i - sin(phi_i)));

ybar = (ybar_o*Ao - ybar_i*Ai)/(Ao-Ai);   % neutral bending plane

%% calculate configuration

h = L*(ro + ybar)/(R * n);   % calc notch height
u = (L - L * (1/R) * (ro + ybar))/n; % calc uncut height







% h = thetaMax * (ro+ybar) / n;
% u = L/n - h;


%% Plot
% theta_max = L/R;
% theta = linspace(0, theta_max, 100);
% 
% x = R - R*cos(theta);
% y = R*sin(theta);
% 
% label = sprintf("notch height h= %.2fmm \nuncut height u= %.2fmm",h,u);
% 
% figure;
% plot(x,y, 'r', 'LineWidth',2)
% axis equal, grid on
% xlabel('X [mm]'), ylabel('Z [mm]');
% ylim([0 L])
% text(R-10,2,label);
% title('Wrist of Given Considerations');

end


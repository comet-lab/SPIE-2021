% Script to calculate the notch geometry
clear, close all, clc

% Define the total notch + uncut section length:
L = 5e-3; %[m]

% Define the desired maximum bending angle for the wrist
rMin = 3 * 10^-3; %[m]
%thetaMax = 1; % [radians]
thetaMax = L/rMin; % [radians]
% Define the number of notches
n = 9;

% Tube diameters
OD = 1.1 * 10^-3; % [m] tube outer diameter
ID = 0.9 * 10^-3; % [m] tube inner diameter
w = 0.85*OD; % [m]

%% You should not need to change anything from the code below

ro = OD/2;         % [m] tube outer radius
ri = ID/2;         % [m] tube inner radius

t = w - ro;   % [m]
phio = 2 * acos(t / ro); % [rad]
phii = 2 * acos(t / ri); % [rad]
ybaro = (4 * ro * (sin(0.5 * phio)) ^ 3)/ (3 * (phio - sin(phio)));
ybari = (4 * ri * (sin(0.5 * phii)) ^ 3)/ (3 * (phii - sin(phii)));
Ao = ( (ro ^ 2) * ( phio - sin(phio))) / 2;
Ai = ( (ri ^ 2) * ( phii - sin(phii))) / 2;
ybar1 = (ybaro * Ao - ybari * Ai) / (Ao - Ai);

fy = @(r,theta) (r.*cos(theta));
fA = @(r,theta) (1);
phi = @(r) acos(t./r);
A = integral2( @(r,theta) r.*fA(r,theta), ri,ro,@(c)-phi(c),@(c)phi(c),'AbsTol',1e-12,'RelTol',1e-12 );
ybar = 1/A*integral2( @(r,theta) r.*fy(r,theta),ri,ro,@(c)-phi(c),@(c)phi(c),'AbsTol',1e-12,'RelTol',1e-9 );



h = thetaMax * (ro+ybar) / n;
%u = (L - L * (1/rMin) * (ro + ybar))/n;
u = L/n - h;
%d = h/(h+u);

disp(['Notch Width (w): ' num2str(w*1e3) ' mm']);
disp(['Notch Height (h) : ' num2str(h*1e3) ' mm']);
disp(['Notch Spacing (u): ' num2str(u*1e3) ' mm']);


%return

% Calculate the notch height and the length of the uncut section u:
% h = thetaMax * (ro + ybar);
% u = (s - h);
%

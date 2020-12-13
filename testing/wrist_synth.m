% clear, close all, clc
% % Define the total notch + uncut section length:
% S = 39.67e-3; %[m]
% % Define the desired maximum bending angle for the notch:
% thetaMax = degtorad(90); % [radians]
% % Define the number of notches
% n = 5;
% % Tube diameters
% OD = 1.60 * 10^-3; % [m] tube outer diameter
% ID = 1.40 * 10^-3; % [m] tube inner diameter
% ro = OD/2;         % [m] tube outer radius
% ri = ID/2;         % [m] tube inner radius
% w = 1.40 * 10^-3; % [m]
% t = w - ro;   % [m]
% phio = 2 * acos(t / ro); % [rad]
% phii = 2 * acos(t / ri); % [rad]
% ybaro = (4 * ro * (sin(0.5 * phio)) ^ 3)/ (3 * (phio - sin(phio)));
% ybari = (4 * ri * (sin(0.5 * phii)) ^ 3)/ (3 * (phii - sin(phii)));
% Ao = ( (ro ^ 2) * ( phio - sin(phio))) / 2;
% Ai = ( (ri ^ 2) * ( phii - sin(phii))) / 2;
% ybar1 = (ybaro * Ao - ybari * Ai) / (Ao - Ai);
% fy = @(r,theta) (r.*cos(theta)); 
% fA = @(r,theta) (1);
% phi = @(r) acos(t./r);
% A = integral2( @(r,theta) r.*fA(r,theta), ri,ro,@(c)-phi(c),@(c)phi(c),'AbsTol',1e-12,'RelTol',1e-12 );
% ybar = 1/A*integral2( @(r,theta) r.*fy(r,theta),ri,ro,@(c)-phi(c),@(c)phi(c),'AbsTol',1e-12,'RelTol',1e-9 );
% %h = sqrt(thetaMax * S * d * (ro+ybar))
% %u = h / d - h
% %n = (h + u) / S
% h = thetaMax * (ro+ybar) / n
% u = S/n - h
% d = h/(h+u)
% % h*1e3, u*1e3
% % return
% % % Calculate the notch height and the length of the uncut section u:
% % h = thetaMax * (ro + ybar);
% % u = (s - h);

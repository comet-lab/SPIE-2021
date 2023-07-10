function L = calc_L(view_ang, R, d_cam, z, h)
% alc_arc_params(view_ang, R, d_cam, z, h) returns length [mm] of a notched joint
% - asssumes that the notched joint is to the left of the camera
%
% view_ang [rad] = angle of camera's field of view (phi in drawings)
% R [m] = radius of the wrist's bending motion
% d_cam [m] = x-dir distance from the center of the wrist to camera 
% z [m] = distance from end of endoscope to bottom of the 1st notch
% h [m] = distance from end of endoscope to camera's optical center
%
% Author: I. Chan <iachan@wpi.edu>
% 
% Last Revision: 5/25/2020 

close all
%% far view angle boundary line equation
m = tan(abs(0.5*(pi - view_ang))); % slope 
x = linspace(0, 2*R, 100);
b = -m * d_cam + h;
y = m*x + b;

%% calculate L
[x1, z1] = linecirc(m, b, R, z, R); % find intersection 
ang_R = atan2(z1(1) - z, x1(1) - R); % ang btwn center of circle and intersection
theta_max = pi - ang_R;
L = R * theta_max; 

%% plot circle and far boundary line (for testing purposes)
% % arc representing the wrist curve
% theta = linspace(ang_R, pi, 100);
% x_cir = R*cos(theta) + R;
% y_cir = R*sin(theta) + z;
% 
% % plot
% figure(1)
% hold on
% plot(x, y, '--', 'LineWidth', 2)
% plot(x_cir, y_cir, 'LineWidth', 2)
% plot(x1(1), z1(1), 'm.', 'MarkerSize', 30)
% legend('View Boundary', 'Notched Joint', 'Intersection', 'Location','southeast')
% axis equal, grid on 
% ylim([0, 50])
% title('Notched Joint Curve constrained by Camera View Angle')
% xlabel('X [mm]'), ylabel('Z [mm]')
% hold off

end
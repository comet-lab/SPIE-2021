function [h, u, L] = calc_all(view_ang, R, n, d_cam, z_cam, h_cam, ro, ri, w)
%CALC_ALL calculates the arc legnth and the wrist configurations given the
%other parameters. Uses calc_L to get L, then uses that as the input to
%calc_config
%
% view_ang [rad] = angle of camera's field of view (phi in drawings)
% R [mm] = radius of the wrist's bending motion
% n = number of notches
% d_cam [mm] = x-dir distance from the center of the wrist to camera 
% z_cam [mm] = distance from end of endoscope to bottom of the 1st notch
% h_cam [mm] = distance from end of endoscope to camera's optical center
% *optional* 
% ro [mm] = outer radius of tube
% ri [mm] = inner radius of tube
% w [mm] = depth of notches
%
% Last Revision: 7/16/2020

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

L = calc_L(view_ang, R, d_cam, z_cam, h_cam);
[h, u] = calc_config(L, R, n, ro, ri, w);

end


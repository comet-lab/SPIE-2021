function config_list = gen_config_list(view_ang, R, d_cam, z, h)
%GEN_CONFIG_LIST generates list of notch configurations
%   view_ang [rad] = angle of camera's field of view (phi in drawings)
%   R [mm] = radius of the wrist's bending motion
%   d_cam [mm] = x-dir distance from the center of the wrist to camera 
%   z [mm] = distance from end of endoscope to bottom of the 1st notch
%   h [mm] = distance from end of endoscope to camera's optical center
%
%   config_list = [L #Notches h u]
% 
% Author: I. Chan <iachan@wpi.edu>
%
% Last Revision: 6/7/2020

% initialize lists
L_list = [];
u_list = [];
h_list = [];
notches_list =[];

% calculate the maximum arc length and bending angle
L_max = calc_L(view_ang, R, d_cam, z, h);
% tolerance
L = [L_max; L_max*0.97; L_max*0.95];
% calculate maximum number of notches
% roughly based h and u of initial tests (h + u = 3mm)
max_notches = floor(min(L) / 3); 

% generate configurations a range of notches for max L
% - Fig.5 of the ISMR Paper (Computational Optimization...)
notches = linspace(5, max_notches, max_notches - 4).';  

for ind = 1:length(L)
    for n = 1: length(notches)
        [h,u] = calc_config(L(ind), R, notches(n));
        h_list = [h_list; h.'];
        u_list = [u_list; u.'];
        L_list = [L_list; L(ind)];
        notches_list = [notches_list; notches(n)];
    end
end

% [L #Notches h u]  
config_list = [L_list notches_list h_list u_list];

end


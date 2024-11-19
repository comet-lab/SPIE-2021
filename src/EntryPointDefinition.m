clear all, clc; 


addpath('kinematics')
addpath('utils')
addpath('utils/stlTools')
addpath('path-planning')
addpath('../anatomical-models')

simulationID = 'larynx1-nowrist-dq-0.06-20pts';
load([simulationID '.mat']);


% Read the configuration file to extract information about the meshes
fid = fopen(fullfile('..', 'anatomical-models', 'configurations.txt'));
text = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

configurations = cell2mat(text(2:end));
line_no = find(strcmp(text{1}, modelID));

path = fullfile('..', 'anatomical-models', modelID);

init_config  = configurations(line_no, 1:6);
entry_point  = configurations(line_no, 7:9);
tip_base     = configurations(line_no, 10:12);

% Calculate the transformation to the space of the CT scan
newZ = tip_base .* 1e-3 - entry_point .* 1e-3;
newZ = newZ ./ norm(newZ);
v = cross([0 0 1], newZ);
R = eye(3) + skew(v) + skew(v)^2 * (1-dot([0 0 1], newZ))/norm(v)^2;
R = R * [0 -1 0; 1 0 0; 0 0 1];
t = entry_point .* 1e-3;
T = [R t'; 0 0 0 1];

% now slide the endoscope back by its length, so that all the different
% designs start exploring from the same point
wrist_len_offset = sum(cutouts.u) + sum(cutouts.h); % len of wrist bend section
endo_len_offset = 13.4e-3 + 28.2e-3;                % tip len + bend section len

Tz = eye(4);
Tz(3,4) = -wrist_len_offset - endo_len_offset;
Tz(3,4) = - endo_len_offset;
T = T * Tz;

% Read the meshes from file
pathStl = fullfile('..', 'anatomical-models', modelID, 'tissue.stl');
[vertices, faces, ~, ~] = stlRead(pathStl);
earModel.vertices = vertices*1e-3;
earModel.faces = faces;
earModel.baseTransform = T;

numFaces = length(faces);
colorMap = zeros(numFaces, 1);

% color indices for painting the larynx
colors = [0.2422    0.1504    0.6603;   % purple
          0.9769    0.9839    0.0805];   % yellow

figure('Name', simulationID);
stlPlot(earModel.vertices, earModel.faces, '', colorMap, colors);
hold on

ii = 1;



% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset;
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
scatter3(T(1,4), T(2,4), T(3,4), 100, 'k', 'filled');
hold on

robot.fwkine(qList(:,ii), T);                   % generate fwkin
robotPhysicalModel = robot.makePhysicalModel(); % generate meshes

% handler for mesh plots
h2 = surf(robotPhysicalModel.surface.Xw, ...
    robotPhysicalModel.surface.Yw, ...
    robotPhysicalModel.surface.Zw, ...
    'FaceColor','blue');
h3 = surf(robotPhysicalModel.surface.Xe, ...
    robotPhysicalModel.surface.Ye, ...
    robotPhysicalModel.surface.Ze, ...
    'FaceColor','red');
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('3D Model in STL')
set(gca,'FontSize',12);

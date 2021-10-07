clear all; clc;

simulationID = 'larynx1-PointsRemoved-dq-0.06-10000pts';

load([simulationID '.mat']);
fid = fopen(fullfile('..', 'anatomical-models', 'configurations.txt'));
text = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

configurations = cell2mat(text(2:end));
line_no = find(strcmp(text{1}, modelID));

path = fullfile('..', 'anatomical-models', modelID);

% Read the Raw Meshes from file
pathMe = fullfile(path, 'tissue_cropped.stl');
[vertices, faces, ~, ~] = stlRead(pathMe);
numverts = size(vertices, 2);

% Convert the raw meshes into objects that can be passed
% to the `patch' function
meMesh.faces = faces;
meMesh.vertices = vertices .* 1e-3;


if exist('visibleMapTotal', 'var')
    v = logical(visibleMapTotal);
else
    v = zeros(length(faces),1);
end

figure('Name', simulationID);
stlPlot(meMesh.vertices * 1e3, meMesh.faces, 'Visibility', v);
axis equal
view(-45, 35)
% larynx 1 az = -45, el = 35
% larynx 2 az = -135, el = 35
zlim([-59 -18])
% larynx 1 zlim zmin = -59, zmax = -18
% larynx 2 zlim zmin = -60, zmax = -15
set(gca, 'fontsize', 18);
title('Steerable Fiber');

clear all; clc; close all;
% include the output of TestSTLfile here!
modelID = 'larynx1';

file2 = 'tissue_closed.stl';

% Obtain and read the files from the directory
%Stl1 = fullfile('..', '..', 'anatomical-models', modelID, file1);
Stl2 = fullfile('..', '..', 'anatomical-models', modelID, file2);
%[vertices,faces,Name1] = stlRead(Stl1);
[vertices,faces,Name2] = stlRead(Stl2);
% n = 10;
% vertices = rand(n, 3)-0.5; % Generate random points
% tetra = delaunayn(vertices); % Generate delaunay triangulization
% faces = freeBoundary(TriRep(tetra,vertices)); % use free boundary as triangulation
n = 10000;
%testp = 2*rand(n,3)-1; % Generate random testpoints
Xblim = 15;
Xtlim = 45;
Yblim = -50;
Ytlim = -10;
Zblim = -70;
Ztlim = 10;
testpX = (Xtlim-Xblim).*rand(n,1) + Xblim; % Generate random testpoints
testpY = (Ytlim-Yblim).*rand(n,1) + Yblim;
testpZ = (Ztlim-Zblim).*rand(n,1) + Zblim;
testp = [testpX testpY testpZ];
in = intriangulation(vertices,faces,testp);
% Plot results
h = trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3));
set(h,'FaceColor','black','FaceAlpha',1/3,'EdgeColor','none');
hold on;
plot3(testp(:,1),testp(:,2),testp(:,3),'b.');
plot3(testp(in==1,1),testp(in==1,2),testp(in==1,3),'ro');

figure
h1 = trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3));
set(h1,'FaceColor','black','FaceAlpha',1/3,'EdgeColor','none');
hold on;
%plot3(testp(:,1),testp(:,2),testp(:,3),'b.');
plot3(testp(in==1,1),testp(in==1,2),testp(in==1,3),'b.');
plot3(testp(in==1,1),testp(in==1,2),testp(in==1,3),'ro');
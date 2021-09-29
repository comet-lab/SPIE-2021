clear all; clc; close all;
% include the output of TestSTLfile here!
modelID = 'larynx1';
file2 = 'tissue_closed.stl';

% Change the file accordingly
load('larynx1dq-0.06-5000pts.mat', 'pList')
testp = pList' * 1000;


% Obtain and read the files from the directory
Stl2 = fullfile('..', '..', 'anatomical-models', modelID, file2);
[vertices,faces,Name2] = stlRead(Stl2);
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
plot3(testp(in==1,1),testp(in==1,2),testp(in==1,3),'b.');
plot3(testp(in==1,1),testp(in==1,2),testp(in==1,3),'ro');

pList = [testp(in==1,1),testp(in==1,2),testp(in==1,3)]';
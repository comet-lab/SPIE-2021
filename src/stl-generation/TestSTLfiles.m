% Plot the STL files to manually determine the entry point (X, Y, Z) in
% the Larynx.

% Select the larynx model you want to visualize acording with the folder
% name.

%clear all; close all; clc;

modelID = 'larynx7a';

%%
% Files that we are going to read
file1 = 'tissue.stl';
file2 = 'tissue_cropped.stl';

% Obtain and read the files from the directory
Stl1 = fullfile('..', '..', 'anatomical-models', modelID, file1);
Stl2 = fullfile('..', '..', 'anatomical-models', modelID, file2);
[V1,F1,Name1] = stlRead(Stl1);
[V2,F2,Name2] = stlRead(Stl2);

%% Plot the STL files 

% Plot the larynx for calcReachableSpace function 
figure
stlPlot(V1,F1,Name1);
title(upper(file1(1:6)));

% Plot the larynx model for calcVisibleArea function
figure
stlPlot(V2,F2,Name2);
title('TISSUE CROPPED');
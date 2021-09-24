clear all; clc;

% Load the file to remove the points 
load('larynx1dq-0.06-5000pts.mat');

otherinfo = [];
otherinfo = [otherinfo '-PointsRemoved-'];
%%
% Generate the new name of the simulation
simulationID = [modelID otherinfo 'dq-' num2str(dq) '-' num2str(nPoints) 'pts'];

%% Set up the threshold limits
% Look into the original plot to determine the threshold limits in the X, Y
% or Z axis as appropriate.

Ly1 = -0.015; % Threshold limits on Y for Larynx 1 (Ly1 and Ly2)
Ly2 = -0.044; 

%% Use of the find function
% Find the position of the points outside the larynx in the pList matrix. 
% The outcome of the find function is the row (R) and the column (C) of
% each reachable point that is outside of the larynx.

%[R,C] = find(pList(2,:) > -0.015); % First test - for one point
%[R,C] = find(pList(2,:) < -0.044); % Second test

[R,C] = find(pList(2,:) > Ly1 | pList(2,:) < Ly2 ); % pList(a,:) where a = 1 for x, 2 for y and 3 for z axis

%%
% Generate the coordinates (x,y,z) matrix of points that are going to
% removed.

pList2remove = zeros(3,length(C));

for i = 1: length(C)
    
    pList2remove(:,i) = pList(:,C(1,i));

end
%% 
% Delete the rows and columns of the he points that are outside the larynx
pList(:,C) = [];

%%
% Save and plot the result
save([simulationID '.mat']);
calcVisibleArea(simulationID, 'mcrc', laserOffsetAngle);
makeVisibilityFig(simulationID);
filename = 'testSim.csv';
getSimData(simulationID, filename);
savefig(['figures/' simulationID '.fig']);

%% Generate simulation video
% Optional, it could be suppressed.
animateResults(simulationID);
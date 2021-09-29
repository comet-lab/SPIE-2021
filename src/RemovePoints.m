clear all; clc;

% Load the file to remove the points 
load('larynx1dq-0.06-5000pts.mat');

otherinfo = [];
otherinfo = [otherinfo '-PointsRemoved-'];
%%
% Generate the new name of the simulation
simulationID = [modelID otherinfo 'dq-' num2str(dq) '-' num2str(nPoints) 'pts'];

% % %% Set up the threshold limits
% % % Look into the original plot to determine the threshold limits in the X, Y
% % % or Z axis as appropriate.
% % 
% % Ly1 = -0.015; % Threshold limits on Y for Larynx 1 (Ly1 and Ly2)
% % Ly2 = -0.044; 
% % 
% % %% Use of the find function
% % % Find the position of the points outside the larynx in the pList matrix. 
% % % The outcome of the find function is the row (R) and the column (C) of
% % % each reachable point that is outside of the larynx.
% % 
% % %[R,C] = find(pList(2,:) > -0.015); % First test - for one point
% % %[R,C] = find(pList(2,:) < -0.044); % Second test
% % 
% % [R,C] = find(pList(2,:) > Ly1 | pList(2,:) < Ly2 ); % pList(a,:) where a = 1 for x, 2 for y and 3 for z axis
% % 
% % %%
% % % Generate the coordinates (x,y,z) matrix of points that are going to
% % % removed.
% % 
% % pList2remove = zeros(3,length(C));
% % 
% % for i = 1: length(C)
% %     
% %     pList2remove(:,i) = pList(:,C(1,i));
% % 
% % end
% % %% 
% % % Delete the rows and columns of the he points that are outside the larynx
% % pList(:,C) = [];

file2 = 'tissue_closed.stl';
testp = pList' * 1000;
Stl2 = fullfile('..', 'anatomical-models', modelID, file2);
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

pList = ([testp(in==1,1),testp(in==1,2),testp(in==1,3)]')/1000;



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
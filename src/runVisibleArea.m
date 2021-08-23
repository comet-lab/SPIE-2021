clc;

% %% Load the mat file that you want to run
% 
% matfile = 'larynx8a-nowrist-dq-0.06-10000pts.mat';
% load(matfile);

%%
% First, load the mat file

%% offste angle [0 45 70 90] deg
useWrist = false;
laserOffsetAngle = 90;

%% Anatomical model definition
modelID = 'larynx8a'; % ID of the anatomical model (see the `anatomical-models' folder)

otherinfo = [];

if ~useWrist
    otherinfo = [otherinfo '-nowrist-'];
end

if laserOffsetAngle
    otherinfo = [otherinfo '-Laser_ang-' num2str(laserOffsetAngle) '-'];
end

simulationID = [modelID otherinfo 'dq-' num2str(dq) '-' num2str(nPoints) 'pts'];

if laserOffsetAngle
    fprintf("Laser Angle offset of %d deg\n", laserOffsetAngle)
end

save([simulationID '.mat']);

%% Run Ray casting
calcVisibleArea(simulationID, 'mcrc', laserOffsetAngle);
%% Create a video of this simulation and Histogram
animateResults(simulationID);
makeVisibilityFig(simulationID);

filename = 'testSim.csv';
getSimData(simulationID, filename);

% Save figure to folder
savefig(['figures/' simulationID '.fig']);
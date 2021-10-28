clear; clc; close all;

%% Load the mat file that you want to run
% First, load the mat file
matfile = 'larynx6b-nowrist-dq-0.06-1000pts.mat';
load(matfile);

%% offste angle 
% select the angle you want to run --> [0 45 70 90] deg
useWrist = false;
laserOffsetAngle = 90;

%% Anatomical model definition

modelID = 'larynx6b'; % ID of the anatomical model (see the `anatomical-models' folder)

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

save([simulationID '.mat']); % save the new workspace with the new offset 
                             % angle and generated the ID for the new 
                             % simulations

%% Run Ray casting
calcVisibleArea(simulationID, 'mcrc', laserOffsetAngle);
%% Create a video of this simulation and Histogram

makeVisibilityFig(simulationID);

filename = 'testSim.csv';
getSimData(simulationID, filename);

% Save figure to folder
savefig(['figures/' simulationID '.fig']);
animateResults(simulationID);
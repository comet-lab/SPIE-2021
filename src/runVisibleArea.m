clear all; clc; close all;

%% Load the mat file that you want to run

load('larynx7a-nowrist-PointsRemoved-dq-0.06-10000pts.mat');

%%
% First, load the mat file

%% offste angle [0 45 70 90] deg
useWrist = false;
laserOffsetAngle = 45;

%% Anatomical model definition

otherinfo = [];

if ~useWrist
    otherinfo = [otherinfo '-nowrist' '-PointsRemoved-'];
end

if laserOffsetAngle
    otherinfo = [otherinfo 'Laser_ang-' num2str(laserOffsetAngle) '-'];
end

simulationID = [modelID otherinfo 'dq-' num2str(dq) '-' num2str(nPoints) 'pts'];

if laserOffsetAngle
    fprintf("Laser Angle offset of %d deg\n", laserOffsetAngle)
end

save([simulationID '.mat']);

%% Run Ray casting
tstart = tic;
calcVisibleArea(simulationID, 'mcrc', laserOffsetAngle);
telapsed = toc(tstart);
s = seconds(telapsed);
s.Format = 'hh:mm:ss'
%% Create a video of this simulation and Histogram
makeVisibilityFig(simulationID);
filename = 'testSim.csv';
getSimData(simulationID, filename);
% Save figure to folder
savefig(['figures/' simulationID '.fig']);
animateResults(simulationID);
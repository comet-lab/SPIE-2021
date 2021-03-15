%% script to make video demonstrating how the gap is produced
close all; clc; clear all;

%% make endoscope and laser 
OD = 0.60e-3;                     % endoscope outer diameter [m]
ID = 0.20e-3; 

n =10;
cutouts.w = 1.36 * ones(1,n) * 1e-3; % [m]
cutouts.u = [0.5 * ones(1,n-1) * 1e-3, (4.5+0.92) * 1e-3]; % [m]
cutouts.h = 0.17 * ones(1,n) * 1e-3; % [m]
cutouts.alpha = zeros(1,n);

robot = EndoWrist(ID, OD, n, cutouts);
robot.endo.bend_sec = 5e-3;
wristConfig = [0, 0, 0];

%% Initial Config
%    k  theta  dz  dl phi tau
q = [0 0 0 wristConfig];
robot.fwkine(q);
robotPhysicalModel = robot.makePhysicalModel();
pts = zeros(3,1);

pts(:,1) = robot.pose();

fig = figure;
hold on; axis equal; grid on;
h1 = scatter3(pts(1,:), pts(2,:), pts(3,:), 'filled');

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

ylim([-.008 .008]);
xlim([-.01 .005]);

view([-90 -75]);
title('Visualization of Gap in Reachable Points')
set(gca,'FontSize',12);
camlight('headlight');
material('dull');

%% iterate through config
maxBend = 0.025; % [1/m]
maxKappa= 1/maxBend;
kinc = 3;
kappas = [0:kinc:40 40:-kinc:0 0:-kinc:-40 -40:kinc:0];
maxTheta = 90;

tsteps = 5;
thetas = [linspace(0, maxTheta, tsteps), linspace(-22.5, -maxTheta, tsteps)];

i = 2;
pauseTime = .1;
video = VideoWriter('gapvisualization_front', 'MPEG-4'); %create the video object
video.FrameRate = 40;
video.Quality = 100;

open(video); %open the file for writing
for t = thetas
    for k = kappas
        q = [k deg2rad(t) 0 wristConfig];
        robot.fwkine(q);
        m = robot.makePhysicalModel();
        pts(:,i) = robot.pose();

        % update scatter plot
        h1.XData = pts(1,:);
        h1.YData = pts(2,:);
        h1.ZData = pts(3,:);

        h2.XData = m.surface.Xw;
        h2.YData = m.surface.Yw;
        h2.ZData = m.surface.Zw;

        h3.XData = m.surface.Xe;
        h3.YData = m.surface.Ye;
        h3.ZData = m.surface.Ze;

        i = i + 1;
        f = gcf;
        f.Position = [100 100 1080 960];    
        frame = getframe(fig);
        
        writeVideo(video, frame);
    end 
    pause(pauseTime*2)
end
close(video);








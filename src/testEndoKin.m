%% Script to test the kinematics of the the Endoscope
clc, clear all
addpath('kinematics')
addpath('path-planning')
addpath('utils')

%% Create the Endoscope

% initial configuration
bend_angle = deg2rad(45);
s = 50e-3;
kappa = bend_angle/s;
theta = 0;
dz = 10e-3;
config = [kappa s theta dz];
% 
% endo = Endoscope();     % create endoscope
% endo.fwkine(config);    % fwkine for initial configuration
% P = endo.pose(:,end);
% T = endo.transformations;
% 
% Pcam = endo.camT(1:3,4);        % pos of camera
% Pwrist = endo.wristT(1:3,4);    % pose of wrist
% 
% endoModel = endo.makePhysicalModel();

%% Create the wrist
n = 7; % number of cutouts

cutouts.w = 1.36 * ones(1,n) * 1e-3; % [m]
cutouts.u = [0.5 * ones(1,n-1) * 1e-3, (4.5+0.92) * 1e-3]; % [m]
cutouts.h = 0.17 * ones(1,n) * 1e-3; % [m]
cutouts.alpha = zeros(1,n);

configuration = [sum(cutouts.h), 0, 0];
% wrist = Wrist(1.4e-3, 1.6e-3, n, cutouts);
% wrist.fwkine(configuration, endo.wristT);
% wristModel = wrist.makePhysicalModel();

%% Create EndoWrist
T = [0.6419   -0.2711    0.7173    0.1116;
   -0.2711    0.7949    0.5429    0.0614;
   -0.7173   -0.5429    0.4367    0.0235;
         0         0         0    1.0000];
         
endowrist = EndoWrist(1.4e-3, 1.6e-3, n, cutouts);
q = [40 20e-3 0 10e-3 configuration];
endowrist.fwkine(q);
m = endowrist.makePhysicalModel();

%% Plot curve
% 
% % tip point
% X = P(1,:);
% Y = P(2,:);
% Z = P(3,:);
% 
% figure('Name', "Endoscope Backbone Curve")
% scatter3(X, Y, Z, 100, 'k', 'filled');  % Tip points
% hold on, axis equal
% 
% % points of curve
% plot3(endoModel.backbone(1,:), ...
%       endoModel.backbone(2,:), ...
%       endoModel.backbone(3,:), ...
%       'LineWidth', 2.5);
%   
% xlabel('X[mm]')
% ylabel('Y[mm]')
% zlabel('Z[mm]')
% title('Backbone curve of endoscope');
% 
% %% Plot mesh endoscope
% figure('Name', 'Endoscope Mesh')
% scatter3(X, Y, Z, 100, 'k', 'filled');
% hold on, axis equal
% 
% % plot mesh surfaces
% X = endoModel.surface.X;
% Y = endoModel.surface.Y;
% Z = endoModel.surface.Z;
X = m.surface.Xe;
Y = m.surface.Ye;
Z = m.surface.Ze;
surf(X, Y, Z, 'FaceColor','red');
hold on, axis equal
% 
% scatter3(Pcam(1), Pcam(2), Pcam(3), 100, 'g', 'filled');        % camera point
% scatter3(Pwrist(1), Pwrist(2), Pwrist(3), 100, 'b', 'filled');  % wrist point
% 
% xlabel('X[mm]')
% ylabel('Y[mm]')
% zlabel('Z[mm]')
% title('Endoscope Mesh with wrist')
% 
% %% Plot Wrist
% P = wrist.pose(:,end);
% 
% X = P(1,:);
% Y = P(2,:);
% Z = P(3,:);
% 
% scatter3(X, Y, Z, 100, 'c', 'filled');  % wrist tip point
% 
% % wrist mesh
X = m.surface.Xw;
Y = m.surface.Yw;
Z = m.surface.Zw;
% 
% X = wristModel.surface.X;
% Y = wristModel.surface.Y;
% Z = wristModel.surface.Z;
surf(X, Y, Z, 'FaceColor','blue');
% 
% %% Plot Camera FoV cone
% 
% % surface points
% X = endoModel.cam.X;
% Y = endoModel.cam.Y;
% Z = endoModel.cam.Z;
% surf(X, Y, Z, 'FaceColor','green', 'FaceAlpha', 0.2);   % camera fov mesh with transparency
% 
% 
% %% Check if tip is within bounds
% endo.tipInBounds(wrist.pose(:,end))



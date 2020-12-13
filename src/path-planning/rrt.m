function [qListNormalized,qList,pList,aList] = rrt(robot, qbounds, earModel, ossiclesModel, nPoints)
% RRT implements the basic Rapidly-Exploring Random Trees algorithm for a
% generic continuum robot
%
% Author: L. Fichera <lfichera@wpi.edu>
%
% Last revision: 3/6/2019

if nargin < 3
    collisionDetection = false;
    nPoints = 1000;
else
    collisionDetection = true;
end

% algorithm parameters
deltaQ = [0.05 0.001 0.01]; % step
maxDispl = qbounds(2);
minDispl = qbounds(1);
maxRot   = qbounds(4);
minRot   = qbounds(3);
maxAdv   = qbounds(6);
minAdv   = qbounds(5);

% initialize the tree and the starting point
qListNormalized = zeros(3, nPoints);
qList = zeros(3, nPoints);
pList = zeros(3, nPoints);
aList = zeros(3, nPoints);


% initialize the base transform
if collisionDetection
    T_robot_in_env = earModel.baseTransform;
else
    T_robot_in_env = eye(4);
end

robot.fwkine(qList(:,1), T_robot_in_env);
T = robot.transformations(:,:,end);
pList(:,1) = T(1:3,4,end);
aList(:,1) = T(1:3,3,end);
% iteratively build the tree
hw = waitbar(0, 'Sampling the configuration space. Please wait...');

jj = 1;

%for ii = 1 : nPoints
while true
    % Generate a random point, identify the closest point in the tree
    % and move towards the new point
    qRand = rand(3,1);
    qNearest = nearestVertex(qRand, qListNormalized, jj);
    qNew = move(qNearest, qRand, deltaQ);
    
    % Scale up the point
    displ = qNew(1) * maxDispl;
    rot   = qNew(2) * maxRot;
    adv   = qNew(3) * (maxAdv - minAdv) + minAdv;
    
    % Check for potential collisions
    robot.fwkine([displ, rot, adv], T_robot_in_env);
    T = robot.transformations(:,:,end);
    
    if collisionDetection
        robotPM = robot.makePhysicalModel();
        
        testpts = [robotPM.surface.X(:) robotPM.surface.Y(:) robotPM.surface.Z(:);
                   robot.pose(:,end)'];
        
        collisionMe = intriangulation(earModel.vertices, ...
            earModel.faces, testpts);
        
        collisionOs = intriangulation(ossiclesModel.vertices, ...
            ossiclesModel.faces, testpts);
        
        collisionMe = sum(collisionMe);
        collisionOs = sum(collisionOs);
        
        if collisionMe > 0 || collisionOs > 0
            disp(['Collision detected. Total number of configurations currently in the tree: ' num2str(jj)]);
            continue;
        end
    end
    
    % If no collision, add this point to the tree
    qListNormalized(:,jj) = qNew;
    
    qList(1,jj) = displ;
    qList(2,jj) = rot;
    qList(3,jj) = adv;
    
    pList(:,jj) = T(1:3,4,end);
    aList(:,jj) = T(1:3,3,end);
    
    jj = jj + 1;
    
    if jj > nPoints
        break;
    end
    
    waitbar(jj/nPoints, hw, 'Sampling the configuration space. Please wait...');
end

close(hw);

qListNormalized = qListNormalized(:,1:jj-1);
qList = qList(:,1:jj-1);
pList = pList(:,1:jj-1);
aList = aList(:,1:jj-1);
end
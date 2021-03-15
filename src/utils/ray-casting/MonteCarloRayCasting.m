function [visibleFaces visible2seen, quiver] = MonteCarloRayCasting(viewPoint, approach, vertices,...
    seenFaces, rangeFaces, nRays, FOV, laserRange)
%MONTECARLORAYCASTING Creates a list of visible faces
% Output:
%   visibleFaces - [3xL] list of faces visible from point
%   seen2visible - [1xL] array of seen index of each visible face
%   quiver - struct containing ray info to be plotted
%       vp - [3xnRays] viewpoint repeated for each ray
%       dispR - [3xnRays] ray coming out of viewpoint scaled to hit the
%       face
%
% Inputs:
%   viewPoint - [3x1] point to look out from
%   approach  - [3x1] direction of viewPoint
%   vertices  - [3xN] list of 3D coords for vertices of mesh
%   seenFaces - [3xM] list of faces that are possibly seen by the tip (rays
%       generated to their centroids)
%   rangeFaces -[3xL] list of all faces within a range of the tip
%   nRays     - number of rays to generate
%   FOV    - [rad] angle of laser
%   laserRange - [mm]
%
%   Notes:
%       Select which version to use by changing the values of 'toFace'
%           toFace will generate rays to the centroid of each face instead
%           of randomly within a cone
%       display will show the result of the rays with the seenFaces
%
%
%   Authors: Jesse F. d'Almeida  <jfdalmeida@wpi.edu>
%
%   Last Revision: 7/24/2020

%% Debug parameters (keep false for normal use)
toFace = false;
display = false;

%% Set up
% scale up size of vertices and viewpoint for visibility
scalar = 1e3;
vertices = vertices * scalar;
viewPoint = viewPoint  * scalar;

% vertices that make up faces to intercept
v1 = vertices(:,rangeFaces(1,:));
v2 = vertices(:,rangeFaces(2,:));
v3 = vertices(:,rangeFaces(3,:));

% vertices that make up faces to generate rays to
c1 = [];
c2 = [];
c3 = [];

seenFaceSize = size(seenFaces,2);
visibleFaces = []; % faces visible on surface
visible2seen = []; % indices of seen faces that are visible
quiv = [];         % holds quiver data for plotting rays

%% Create Rays
% toFace is true for a ray to the avg of each face, false to randomly
% generate rays
if toFace
    rays = zeros(3,seenFaceSize);
    len = seenFaceSize;
    
    % vertices that make up faces to generate rays to
    c1 = vertices(:,seenFaces(1,:));
    c2 = vertices(:,seenFaces(2,:));
    c3 = vertices(:,seenFaces(3,:));
else
    rays = RandomRays(rad2deg(FOV/2), approach, nRays);
    rays = normc(rays);
    len = nRays;
end
rayScalar = zeros(1,len);

%% Find closest visible face
for i = 1:len
    if toFace
        % make ray from viewpoint to the centroid of face
        vs = [c1(:,i)-viewPoint, c2(:,i)-viewPoint, c3(:,i)-viewPoint];
        r = mean(vs, 2);
        rays(:,i) = r;
    end
    ray = rays(:,i);
    
    % use triangle/ray intersect to get dists of intersection between ray
    % and faces
    [inter,dists, ~, ~, ~] = TriangleRayIntersection(...
        viewPoint', ray', v1, v2, v3,...
        'border', 'normal');
    
    % if there is an intersection
    flagged = sum(inter) > 0;
    if flagged
        dists(dists <=0) = inf;
        [val1, idx1] = min(dists);
        
        dists(inter == 0)= inf;  % only want the distances where intercepts exist
        [val2, idx2] = min(dists);  % min dist is visible face to the ray
        
        val = val2;
        idx = idx2;
        
        % only add if within laser range
        if val <= laserRange
            visibleFaces(:,end+1) = rangeFaces(:, idx);    % add face at index
            visible2seen(end+1) = idx;      % add index of the seen face
            rayScalar(i) = val;             % store dist of interception
        end
    end
end

%% Display Rays (can comment out)
% displays the seen faces of the mesh with the rays pointing from viewpoint
rayScalar(rayScalar == 0) = max(rayScalar);
dispR = rays .* rayScalar;  % scale ray lengths by their intersection distance
l = size(dispR,2);
vp = repmat(viewPoint,1,l);
quiv = [vp;dispR];
quiver.vp = vp;
quiver.dispR = dispR;
if display
    figure
    subplot(121)
%     stlPlot(vertices' , seenFaces', '');
    title('Seen faces with rays')
%     hold on, axis equal
    quiver3(vp(1,:), vp(2,:), vp(3,:), rays(1,:), rays(2,:), rays(3,:));
    
    subplot(122)
    stlPlot(vertices' , rangeFaces', '');
    hold on, axis equal
    title('Faces within range with rays')
    quiver3(vp(1,:), vp(2,:), vp(3,:), dispR(1,:), dispR(2,:), dispR(3,:));
end
end


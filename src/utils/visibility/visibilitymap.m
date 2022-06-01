function [visibleMap, quiver] = visibilitymap(viewPoint, approachVec, meModel, alg)
% VISIBILITYMAP - list of visible areas from a specific viewpoint within a
% model and a given approach
%   
% Outputs:
%   visibleMap:   (1xN) boolean array of if each face is visible
%   quiver - struct containing ray info to be plotted
%       vp - [3xnRays] viewpoint repeated for each ray
%       dispR - [3xnRays] ray coming out of viewpoint scaled to hit the
%       face
%
% Inputs:
%   viewPoint: (3x1) point to look out from
%   approachVec: (3x1) direction for tip form viewpoint
%   meModel: (model obj) primary model 
%   alg: ('hpr' or 'mcrc') to select between which algorithm is used (Hidden
%       Point Removal or Monte Carlo Ray Casting
%
%   Authors: Jesse F. d'Almeida  <jfdalmeida@wpi.edu>
%   
%   Last Revision: 6/16/2020

% Laser Parameters
laserRange = 10e-3; %Return to 3e-3;      % range of laser is 3mm
%laserFOV = 40;          % divergence angle of laser Endostat Fiber
laserFOV = 85; %return to 60;          % divergence angle of laser Optical Fiber FP200ERT Thorlabs
FOV = deg2rad(laserFOV);

% model properties
faces = [meModel.faces]';
vertices = [meModel.vertices]';

%% set outputs
visibleMap = zeros(size(faces,2),1);
quiver.vp = [];
quiver.dispR = [];

%% Find seen vertices
% Cast rays from the viewPoint to each of the centroids
rays = bsxfun(@minus, vertices, viewPoint);

% rays of length that are within the laser range
rayLength = vecnorm(rays);
closeVertices = (rayLength <= laserRange);
rangeVertices = (rayLength <= 7e-3);

% Convert rays into unit vectors
raysu = normc(rays);

% Make as many copies of approachVec as the number of rays we generated
approachVecRep = repmat(approachVec, 1, length(raysu));

% See what rays fall within the "field of view" of the camera
product = sum(approachVecRep .* raysu);

% get vertices between the bounds
fovVertices = (product >= cos(FOV/2));

% array [n] where each index is boolean of if vertex is seen
seenVertices = closeVertices .* fovVertices;

% exit is no vertices are seen
if ~logical(seenVertices)
    return
end

%% Get list of seen faces that are made up of the seen vertices
% loop through each face to find which are made up of seen vertices
seenFaces = [];
seen2Faces = [];   % holds original idx of each seen face
anyFacesSeen = false;

for ii = 1:size(faces,2)
    sf = faces(:,ii);
    v1 = sf(1);
    v2 = sf(2);
    v3 = sf(3);
    
    % if each vertex in face is 'seen', add to list
    if seenVertices(v1) && seenVertices(v2) && seenVertices(v3)
        anyFacesSeen = true;
        seenFaces(:,end+1) = sf;     % add face to the end of matrix
        seen2Faces(end+1) = ii;      % index of face that is seen
    end
end


%% Get list of all the faces within a range of the tip
rangeFaces = [];
range2Faces = [];   % holds original idx of each face in range

for ii = 1:size(faces,2)
    sf = faces(:,ii);
    v1 = sf(1);
    v2 = sf(2);
    v3 = sf(3);
    
    % if each vertex in face is within the range, add to list
    if rangeVertices(v1) && rangeVertices(v2) && rangeVertices(v3)
        rangeFaces(:,end+1) = sf;     % add face to the end of matrix
        range2Faces(end+1) = ii;      % index of face that is seen
    end
end

%% Find visible faces ray casting
if strcmp(alg,'mcrc')
    nRays = 1000;
    [~, visible2seen, quiver] = MonteCarloRayCasting(viewPoint, approachVec,...
    vertices, seenFaces, rangeFaces, nRays, FOV, laserRange * 1e3);
    visibleMap(range2Faces(visible2seen)) = 1;

elseif strcmp(alg,'hpr')
%% Find visible faces HPR
    onlySeenVerts = vertices(:, seenVertices == 1);
    visibleVerticesIdx = HPR(onlySeenVerts', viewPoint', 2.5);
    
    seen2Vertices = find(seenVertices);     % idxs of vertices that are seen
    
    lenSeenFaces = size(seenFaces,2);

    % find which of seen faces is visible
    for i = 1:lenSeenFaces
         sf = seenFaces(:,i);

         % verts in the seen face are visible
         if ismember(sf, seen2Vertices(visibleVerticesIdx))
             visibleMap(seen2Faces(i)) = 1;     % convert seen idx to faces and mark as visible
         end
    end
else
    disp('ERROR: not valid algorithm')
end
end

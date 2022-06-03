function calcVisibleArea2(simulationID, alg)
%% Load in Simulation
fprintf('* Estimation of the visible surface. Using %s algorithm*\n', alg)

load([simulationID '.mat']);
fid = fopen(fullfile('..', 'anatomical-models', 'configurations.txt'));
text = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

configurations = cell2mat(text(2:end));
line_no = find(strcmp(text{1}, modelID));

path = fullfile('..', 'anatomical-models', modelID);

% Read the Raw Meshes from file
pathMe = fullfile(path, 'tissue_cropped.stl');
[vertices, faces, ~, ~] = stlRead(pathMe);

% Convert the raw meshes into objects that can be passed
% to the `patch' function
meMesh.faces = faces;
meMesh.vertices = vertices .* 1e-3;
% meMesh.Vertices = meMesh.Vertices ;

meMesh.FaceVertexCData = ones(size(meMesh.vertices, 1), 1);
meMesh.LineStyle = 'none';
meMesh.FaceColor = 'flat';
meMesh.FaceAlpha = 0.4 ;

% Calculate the visual range- list of faces [3xM] that are visible
visibleFaces = [];

% empty quiver for when no faces are visible
q1.vp = [0;0;0];
q1.dispR = [0;0;0];

% empty outputs to fill
quivCamera = repmat(q1, 1, nPoints);          
visibleMapCamera = zeros(length(faces),1);

qList = qList;
T = T;
robot = robot;

hw = waitbar(0, 'Calculating the visibility map. Please wait...');

%% RUN VISIBILITY: for each point, get list of visibile faces

% for each point, calculate visibility 
for jj = 1 : size(pList, 2)

    % First calculate the position of the tip of the endoscope
    %endoTip = robot.endo.transformations(1:3,4,end);
    %endoApproach = robot.endo.transformations(1:3,4,end);
    
    % get current transformation to tip
    robot.fwkine(qList(:,jj), T);
    endoTip = robot.endo.camT(1:3,4,end);
    approach = robot.endo.transformations(1:3,3,end);
    %approach = approach4(1:3);                      % split 3D vector

    [m, q] = visibilitymap(endoTip, approach, meMesh, alg, 50e-3, 85, visibleMapTotal);
    visibleMapCamera(:,jj) = m;   % visibility of faces for this point
    
    % add if there are quivers, if not add empty so no rays will be
    % displayed
    if q.vp
        quivCamera(jj) = q;
    else
        quivCamera(jj) = q1;
    end
    
    waitbar(jj/nPoints, hw, 'Calculating the visibility map. Please wait...');
end
close(hw);

% Sum rows of map to get amount for each face
numFaces = length(faces);
visibleMapTotalCamera = sum(visibleMapCamera, 2);            % sum rows of visible map
visFACESCamera = sum(logical(visibleMapTotalCamera));        % total visible face
percVISFACESCamera = (visFACESCamera/numFaces) * 1e2;

%% GATHER RESULTS: calculate area visible
visible_areaCamera = seenArea(meMesh, visibleMapTotalCamera);
allFaces = ones(numFaces,1);
total_area = seenArea(meMesh, allFaces);
percVISAREACamera = visible_areaCamera/total_area * 100;
visAREACamera = visible_areaCamera * 1e6;

fprintf('Faces Visible: %d \tPercent of Total Faces: %.2f%% \n',...
    visFACESCamera, percVISFACESCamera)
fprintf('Visible Surface Area: %.2f mm^2 \tPercent of Total Surface Area: %.2f%%\n\n',...
    visAREACamera, percVISAREACamera);

%% OUTPUTS
results.visFacesCamera = visFACESCamera;
results.percFacesCamera = percVISFACESCamera;
results.visAreaCamera = visAREACamera;
results.percAreaCamera = percVISAREACamera;

save([simulationID '.mat'], '-v7.3');


end
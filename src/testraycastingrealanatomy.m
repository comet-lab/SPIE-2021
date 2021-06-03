%% Call testpathplanningrealanatomy first!

close all

% Load the segmentation of the ME
load(fullfile(path, 'record.mat'));

% Read the Raw Meshes from file
pathMe = fullfile('anatomical-models', modelID, 'me.mesh');
pathOs = fullfile('anatomical-models', modelID, 'ossicle.mesh');
rawMeMesh = meshread(pathMe);
rawOsMesh = meshread(pathOs);

% Convert the raw meshes into objects that can be passed
% to the `patch' function
meMesh.Faces = rawMeMesh.triangles' + 1;
meMesh.Vertices = bsxfun(@times, rawMeMesh.vertices', voxel_size);
meMesh.Vertices = meMesh.Vertices .* 1e-3;

osMesh.Faces = rawOsMesh.triangles' + 1;
osMesh.Vertices = bsxfun(@times, rawOsMesh.vertices', voxel_size);
osMesh.Vertices = osMesh.Vertices .* 1e-3;

meMesh.FaceVertexCData = ones(size(meMesh.Vertices, 1), 1);
meMesh.LineStyle = 'none';
meMesh.FaceColor = 'flat';
meMesh.FaceAlpha = 0.4 ;

osMesh.FaceVertexCData = ones(size(osMesh.Vertices, 1), 1);
osMesh.FaceColor = 'flat';
osMesh.FaceAlpha = 0.4 ;


% Calculate the visual range
seenMap = false(rawMeMesh.numverts, length(pList));
ii = 1;
seenMap(:,ii) = visibilitymap(pList(:,ii), aList(:,ii), meMesh, osMesh);

figure('units','normalized','outerposition',[0 0 1 1])
%h1 = patch(meMesh);
%stlPlot(earModel.vertices, earModel.faces, 'Ray Casting test.', seenMap);
%hold on
%patch(osMesh);


%robotPhysicalModel = robot.makePhysicalModel();
%h2 = surf(robotPhysicalModel.surface.X, ...
%    robotPhysicalModel.surface.Y, ...
%    robotPhysicalModel.surface.Z, ...
%    'FaceColor','blue');

%axis equal

%while true
parfor ii = 1 : length(pList)
    seenMap(:,ii) = visibilitymap(pList(:,ii), aList(:,ii), meMesh, osMesh);
    ii
    %robot.fwkine(qList(:,ii), T);
    %robotPhysicalModel = robot.makePhysicalModel();
    
    %colorMap = zeros(length(h1.FaceVertexCData), 1);
    %colorMap(logical(seenMap)) = 5;
    %h1.FaceVertexCData = colorMap;
    
    %h2.XData = robotPhysicalModel.surface.X;
    %h2.YData = robotPhysicalModel.surface.Y;
    %h2.ZData = robotPhysicalModel.surface.Z;
    %title(['Pose ' num2str(ii) ' of ' num2str(size(pList, 2))]);
    
    %fprintf(['Seen area: ' num2str(seenArea(earModel, seenMap)) ' m2.\n']);
    %fprintf('Press "n" to move forward or "p" to move back.\n')
    %fprintf('Press any other key to exit.\n\n')
    
%     while ~waitforbuttonpress, end
%     k = get(gcf, 'CurrentCharacter');
%     
%     switch k
%         case 'p'
%             ii = ii - 1;
%             if ii < 1, ii = 1; end
%         case 'n'
%             ii = ii + 1;
%             if ii > size(pList, 2), ii = size(pList, 2); end
%         case '-'
%             ii = ii + 10;
%             if ii < 1, ii = 1; end
%         case '+'
%             ii = ii + 10;
%             if ii > size(pList, 2), ii = size(pList, 2); end
%         otherwise
%             break
%     end
end

seenMap = sum(seenMap, 2);
seenMap(seenMap > 1) = 1;

figure
h1 = patch(meMesh);
hold on, axis equal, grid on
patch(osMesh);

colorMap = zeros(length(h1.FaceVertexCData), 1);
colorMap(logical(seenMap)) = 5;
h1.FaceVertexCData = colorMap;

%visibility = area_stat(recordnew(:), seenMap(:));

fprintf('Testing complete.\n')
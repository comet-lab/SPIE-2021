function makeVisibilityFig(simulationID, options)
% Produces figures were the visible area of a model is highlighted
% 
% simulationID = identifier of this simulation
%
% Author: Jesse F. d'Almeida  <jfdalmeida@wpi.edu>
arguments
    simulationID (1,:) char
    options.plotVisible (1,1) logical = true
    options.plotCones (1,1) logical = false
    options.plotPoints (1,1) logical = false;
    options.tissueName (1,:) char = 'tissue_cropped'
    options.colorMap = [0 1 1; 0 1 0];
end

load([simulationID '.mat']);
fid = fopen(fullfile('..', 'anatomical-models', 'configurations.txt'));
text = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

configurations = cell2mat(text(2:end));
line_no = find(strcmp(text{1}, modelID));

path = fullfile('..', 'anatomical-models', modelID);

% Read the Raw Meshes from file
pathMe = fullfile(path, [options.tissueName '.stl']);
[vertices, faces, ~, ~] = stlRead(pathMe);
numverts = size(vertices, 2);

% Convert the raw meshes into objects that can be passed
% to the `patch' function
meMesh.faces = faces;
meMesh.vertices = vertices .* 1e-3;

if exist('visibleMapTotal', 'var')
    v = logical(visibleMapTotal);
else
    v = zeros(length(faces),1);
end

if options.plotVisible
%     figure('Name', simulationID);
    gca;
    if options.plotPoints
        subplot(121);
        stlPlot(meMesh.vertices * 1e3, meMesh.faces, 'Visibility', v);
    %     hold on, axis equal
    %     scatter3(collLocs(:,1)*1e3, collLocs(:,2)*1e3, collLocs(:,3)*1e3, 'filled', 'green');

        subplot(122);
        stlPlot(meMesh.vertices * 1e3, meMesh.faces, 'Visibility', v);
        hold on, axis equal
        scatter3(pList(1,:)*1e3, pList(2,:)*1e3, pList(3,:)*1e3, 'filled', 'red');
    else
        stlPlot(meMesh.vertices * 1e3, meMesh.faces, 'Visibility', v, options.colorMap);
    end
else
    v = zeros(length(faces),1);
    stlPlot(meMesh.vertices * 1e3, meMesh.faces, 'Visibility', v);
    hold on, axis equal
    scatter3(pList(1,:)*1e3, pList(2,:)*1e3, pList(3,:)*1e3, 'filled', 'red');
end
if options.plotCones
    for i = 1:size(quiv,2)
        vp = quiv(1,i).vp;
        dispR = quiv(1,i).dispR;
        gcf;
        subplot(122)
        hold on, axis equal
        title('Faces within range with rays')
        quiver3(vp(1,:), vp(2,:), vp(3,:), dispR(1,:), dispR(2,:), dispR(3,:));
    end
end

end
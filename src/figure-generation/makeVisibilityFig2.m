function makeVisibilityFig(simulationID, meMesh, visibleMapTotal, visibleMapTotalCamera)
% Produces figures were the visible area of a model is highlighted
%
% simulationID = identifier of this simulation
%
% Author: Loris Fichera <lfichera@wpi.edu>
%         Based on code by Jesse d'Almeida

% load([simulationID '.mat']);

figure('Name', simulationID);
subplot(131);
v = logical(visibleMapTotal);
stlPlot(meMesh.vertices * 1e3, meMesh.faces, 'Visibility', v);
title('Laser Reachability');

subplot(132);
v = logical(visibleMapTotalCamera);
stlPlot(meMesh.vertices * 1e3, meMesh.faces, 'Visibility', v);
title('Reachability & Visibility');

subplot(133);
v = logical(visibleMapTotal) - logical(visibleMapTotalCamera);
stlPlot(meMesh.vertices * 1e3, meMesh.faces, 'Visibility', v);
title('Difference');

end

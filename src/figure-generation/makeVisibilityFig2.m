function makeVisibilityFig2(simulationID, meMesh, visibleMapTotal, visibleMapTotalCamera)
% Produces figures were the visible area of a model is highlighted
%
% simulationID = identifier of this simulation
%
% Author: Loris Fichera <lfichera@wpi.edu>
%         Based on code by Jesse d'Almeida

% load([simulationID '.mat']);

f = figure('Name', simulationID);
f.Position = [1474, 152, 1090, 643];

c = [[133 122 189] ./ 255;
     1 0 0];

subplot(131);
v = logical(visibleMapTotal);
stlPlot(meMesh.vertices * 1e3, meMesh.faces, 'Visibility', v);
title('(a) Laser Fiber Reachability');
set(gca,'FontSize',16);

subplot(132);
v = logical(visibleMapTotalCamera);
stlPlot(meMesh.vertices * 1e3, meMesh.faces, 'Visibility', v);
title({'(b) Laser Fiber Reachability +','Endoscope Visibility'});
set(gca,'FontSize',16);

subplot(133);
v = logical(visibleMapTotal) - logical(visibleMapTotalCamera);
stlPlot(meMesh.vertices * 1e3, meMesh.faces, 'Visibility', v, c);
title('Difference between (a) and (b)');
set(gca,'FontSize',16);

end

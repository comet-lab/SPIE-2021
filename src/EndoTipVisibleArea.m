simulationID = 'larynx1-nowrist-dq-0.06-11pts';
load([simulationID '.mat']);

% Initialize the variables for the endoscope Tip
TListEndo = zeros(4,4,length(qList));
pListEndo = zeros(3,length(qList));
aListEndo = zeros(3,length(qList));

for i = 1%:length(qList)
   robot.endo.fwkine(qList(1:3,i),T);
   TListEndo(:,:,i) = robot.endo.camT;
   pListEndo(:,i) = TListEndo(1:3,4,i);
   aListEndo(:,i) = TListEndo(1:3,3,i);  
end

results = calcVisibleArea(earModel, ...
                          pListEndo, ...
                          aListEndo, ...
                          'mcrc', ...
                          0, ....
                          robot.endo.cam_range, ...
                          robot.endo.fov);
                      
visibleMap = results.visibleMap;
                      
results = calcVisibleArea(earModel, ...
                          pList, ...
                          aList, ...
                          'mcrc', ...
                          0, ....
                          3e-3, ...
                          40 * pi / 180);
                      
reachableMap = results.visibleMap;

save([simulationID '.mat'], '-v7.3');
load([simulationID '.mat']);


%calcVisibleArea(simulationID, 'mcrc', laserOffsetAngle, endo_range, endo_FOV);
makeVisibilityFig(simulationID);
%animateResults(simulationID);
%% testing shape
Stl1 = fullfile('sphere.stl'); %load stl
[V1,F1,Name1] = stlRead(Stl1);
p = [ 0.8       0.8       1.0;    % green
      0.2422    0.1504    0.6603; % purple
      0.9769    0.9839    0.0805]; % yellow
 figure
stlPlot(V1,F1,Name1,p(2,:))

%% intriangulation function
% n = 10;
% vertices = rand(n, 3)-0.5; % Generate random points
% tetra = delaunayn(vertices); % Generate delaunay triangulization
% faces = freeBoundary(TriRep(tetra,vertices)); % use free boundary as triangulation

n = 1000;
% testp = 2*rand(n,3)-1; % Generate random testpoints
testp = 60*rand(n,3)-1; % Generate random testpoints (change to size of x
in = intriangulation(V1,F1,testp);
% Plot results
h = trisurf(F1,V1(:,1),V1(:,2),V1(:,3));
set(h,'FaceColor','black','FaceAlpha',1/3,'EdgeColor','none');
hold on;
plot3(testp(:,1),testp(:,2),testp(:,3),'b.');
plot3(testp(in==1,1),testp(in==1,2),testp(in==1,3),'ro');

%% remove outside points 
% Generate the coordinates (x,y,z) matrix of points that are going to
% removed.
C = find(in == 0);
testp_remove = zeros(length(C),3);

for i = 1:length(C)
    
    testp_remove(i,:) = testp(C(i),:);

end
%% delete rows
% Delete the rows and columns of the points that are outside the larynx
testp(C,:) = [];
%% plot new point set (without points outside)
figure
h = trisurf(F1,V1(:,1),V1(:,2),V1(:,3));
set(h,'FaceColor','black','FaceAlpha',1/3,'EdgeColor','none');
hold on;
plot3(testp(:,1),testp(:,2),testp(:,3),'b.');

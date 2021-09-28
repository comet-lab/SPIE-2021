

% include the output of TestSTLfile here!

% n = 10;
% vertices = rand(n, 3)-0.5; % Generate random points
% tetra = delaunayn(vertices); % Generate delaunay triangulization
% faces = freeBoundary(TriRep(tetra,vertices)); % use free boundary as triangulation
% n = 1000;
testp = 2*rand(n,3)-1; % Generate random testpoints
in = intriangulation(vertices,faces,testp);
% Plot results
h = trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3));
set(h,'FaceColor','black','FaceAlpha',1/3,'EdgeColor','none');
hold on;
plot3(testp(:,1),testp(:,2),testp(:,3),'b.');
plot3(testp(in==1,1),testp(in==1,2),testp(in==1,3),'ro');
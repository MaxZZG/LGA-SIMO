
clc
close all
clear all

n = 4;
m = 4;
CP = zeros(n+1,m+1,4);

CP(1,1,:) = [0,0,1,1];
CP(2,1,:) = [1,0,0,1];
CP(3,1,:) = [2,0,0,1];
CP(4,1,:) = [3,0,0,1];
CP(5,1,:) = [4,0,0,1];

CP(1,2,:) = [0,1,0,1];
CP(2,2,:) = [1,1,1,1];
CP(3,2,:) = [2,1,1,1];
CP(4,2,:) = [3,1,1,1];
CP(5,2,:) = [4,1,0,1];

CP(1,3,:) = [0,2,0,1];
CP(2,3,:) = [1,2,1,1];
CP(3,3,:) = [2,2,1,1];
CP(4,3,:) = [3,2,1,1];
CP(5,3,:) = [4,2,0,1];

CP(1,4,:) = [0,3,0,1];
CP(2,4,:) = [1,3,1,1];
CP(3,4,:) = [2,3,1,1];
CP(4,4,:) = [3,3,1,1];
CP(5,4,:) = [4,3,0,1];

CP(1,5,:) = [0,4,0,1];
CP(2,5,:) = [1,4,0,1];
CP(3,5,:) = [2,4,0,1];
CP(4,5,:) = [3,4,0,1];
CP(5,5,:) = [4,4,0,1];

Tsurface = T_Surface(1,CP);
plotSurf(1,Tsurface);
Tsurface.plotmesh;
% Tsurface.plotmeshTp;
% Tsurface.plotmeshp;

Tsurface1 = T_Surface(0,CP);
Tsurface1.plotSurf;
Tsurface1.plotmesh;


daspect([1 1 1])
axis equal


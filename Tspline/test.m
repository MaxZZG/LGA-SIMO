
clc
clear all
close all

n = 3;
m = 3;
CP = zeros(n+1,m+1,4);

CP(1,1,:) = [0,0,0,1];
CP(2,1,:) = [1,0,0,1];
CP(3,1,:) = [2,0,0,1];
CP(4,1,:) = [3,0,0,1];

CP(1,2,:) = [0,1,0,1];
CP(2,2,:) = [1,1,1,1];
CP(3,2,:) = [2,1,1,1];
CP(4,2,:) = [3,1,0,1];

CP(1,3,:) = [0,2,0,1];
CP(2,3,:) = [1,2,1,1];
CP(3,3,:) = [2,2,1,1];
CP(4,3,:) = [3,2,0,1];

CP(1,4,:) = [0,3,0,1];
CP(2,4,:) = [1,3,0,1];
CP(3,4,:) = [2,3,0,1];
CP(4,4,:) = [3,3,0,1];


Tsurface = T_Surface(CP);

Tsurface.plotSurf;
Tsurface.plotmesh;
Tsurface.plotmeshTp;
plotDomainCP(Tsurface,1);

Tsurface.plotmeshp;


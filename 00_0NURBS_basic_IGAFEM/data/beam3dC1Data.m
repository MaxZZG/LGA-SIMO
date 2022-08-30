% NURBS data for 3D beam
% NURBS data: control points, knots, basis orders and weights
% Vinh Phu Nguyen, Johns Hopkins University

% control points

a = 10;
b = 4;
c = 2;

noPtsX = 21;
noPtsY = 7;
noPtsZ = 5;

[controlPts,elementVV]=makeB8mesh(a,b,c,noPtsX,noPtsY,noPtsZ);

% basis order

p = 2;
q = 2;
r = 2;

% knot vectors

knotUTemp = linspace(0,1,noPtsX-p+1);
knotVTemp = linspace(0,1,noPtsY-q+1);
knotWTemp = linspace(0,1,noPtsZ-r+1);

uKnot = [0 0 knotUTemp 1 1];
vKnot = [0 0 knotVTemp 1 1];
wKnot = [0 0 knotWTemp 1 1];

% weights

weights = ones(1,noPtsX*noPtsY*noPtsZ)';




          

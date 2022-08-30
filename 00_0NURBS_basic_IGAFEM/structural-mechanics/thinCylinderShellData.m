% Data for the thin cylinder shell problem
% Vinh Phu Nguyen, March 2013
% Cardiff University

addpath ../nurbs-geopdes/inst/
addpath ../post-processing/
addpath ../

R      = 1;
L      = 5;

%% control points

controlPts = zeros(4,9,2);

% first zKnot 

controlPts(1:3,1,1) = [R;0;0];
controlPts(1:3,2,1) = [R;R;0];
controlPts(1:3,3,1) = [0;R;0];
controlPts(1:3,4,1) = [-R;R;0];
controlPts(1:3,5,1) = [-R;0;0];
controlPts(1:3,6,1) = [-R;-R;0];
controlPts(1:3,7,1) = [0;-R;0];
controlPts(1:3,8,1) = [R;-R;0];
controlPts(1:3,9,1) = [R;0;0];

% second zknot

controlPts(1:3,1,2) = [R;0;L];
controlPts(1:3,2,2) = [R;R;L];
controlPts(1:3,3,2) = [0;R;L];
controlPts(1:3,4,2) = [-R;R;L];
controlPts(1:3,5,2) = [-R;0;L];
controlPts(1:3,6,2) = [-R;-R;L];
controlPts(1:3,7,2) = [0;-R;L];
controlPts(1:3,8,2) = [R;-R;L];
controlPts(1:3,9,2) = [R;0;L];

controlPts(4,:,:)   = 1;

fac                 = 1/sqrt(2);

controlPts(4,2,1) = fac;
controlPts(4,4,1) = fac;
controlPts(4,6,1) = fac;
controlPts(4,8,1) = fac;
controlPts(4,2,2) = fac;
controlPts(4,4,2) = fac;
controlPts(4,6,2) = fac;
controlPts(4,8,2) = fac;

% homogenous coordinates (x*w,y*w,z*w)

controlPts(1:3,2,1) = controlPts(1:3,2,1)*fac;
controlPts(1:3,4,1) = controlPts(1:3,4,1)*fac;
controlPts(1:3,6,1) = controlPts(1:3,6,1)*fac;
controlPts(1:3,8,1) = controlPts(1:3,8,1)*fac;
controlPts(1:3,2,2) = controlPts(1:3,2,2)*fac;
controlPts(1:3,4,2) = controlPts(1:3,4,2)*fac;
controlPts(1:3,6,2) = controlPts(1:3,6,2)*fac;
controlPts(1:3,8,2) = controlPts(1:3,8,2)*fac;

%% knot vectors 4x1 mesh

uKnot = [0 0 0 1 1 2 2 3 3 4 4 4];
vKnot = [0 0 1 1];

uKnot = uKnot/max(uKnot);

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

%% refinement at two ends
newKnots  = {[] [0.2 0.8]};
solid     = nrbkntins(solid,newKnots);  
uKnot     = cell2mat(solid.knots(1));
vKnot     = cell2mat(solid.knots(2));
    
%% k-refinement

% p-refinement to raise order

solid = nrbdegelev(solid,[1 2]); 

% then h-refinement 

refineCount = 3;

for i=1:refineCount
    uKnotVectorU = unique(uKnot);
    uKnotVectorW = unique(vKnot);
    
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsZ = uKnotVectorW(1:end-1) + 0.5*diff(uKnotVectorW);
    newKnots  = {newKnotsX newKnotsZ};
    
    % h-refinement
    
    solid     = nrbkntins(solid,newKnots);    
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));
end

% crv = nrbcirc(R,[],deg2rad(0),deg2rad(180));
% nrbplot(crv,80);
% crv = nrbcirc(R-t,[],deg2rad(0),deg2rad(180));
% nrbplot(crv,80);

%%%%%%%%
%% convert NURBS data back to our data structure for analysis

convert2DNurbsShell

%% Material properties

E  = 3e6;
nu = 0.3;
t  = 0.01; % thickness

%% Boundary condition

pressure = 10;


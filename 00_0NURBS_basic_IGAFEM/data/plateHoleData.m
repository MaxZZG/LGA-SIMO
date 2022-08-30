% data for plate with a hole Ck elements
% h and k-refinement can be used.
% Used for a Finite Cell Method.
%
% Vinh Phu Nguyen
% Johns Hopkins University

%%

addpath ../fem_util
addpath ../nurbs-geopdes/inst/
addpath ../nurbs-util/
addpath ../meshing/
addpath ../fem-functions/
addpath ../post-processing/
addpath ../xiga/

%clear all

L = 5; % half plate width

p = 1;
q = 1;

noPtsX = 31;
noPtsY = 62;

gcoord=meshRectangularCoord(L,2*L,noPtsX-1,noPtsY-1);
controlPts=gcoord;

weights = ones(1,noPtsX*noPtsY)';

knotUTemp = linspace(0,1,noPtsX);
knotVTemp = linspace(0,1,noPtsY);

uKnot = [0 knotUTemp 1];
vKnot = [0 knotVTemp 1];

noCtrPts       = noPtsX   * noPtsY;
noDofs         = noCtrPts * 2;

%% generate element connectivity ...

generateIGA2DMesh

% plot the mesh

buildVisualizationMesh;


% circular hole data

r  = 1; % radius
xc = 0.5*L; % x coord of center
yc = L; % y coord of center

VOID = [xc yc r];

levelSetVoids

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(node,elementV,'Q4','b-',1.2);
% plot the circle
theta = 0:0.01:2*pi;
xo = xc + r*cos(theta) ;
yo = yc + r*sin(theta) ;
plot(xo,yo,'k-','Linewidth',1.9);
% plot elements cut by the circle
plot_mesh(node,elementV(splitElems(:,1),:),'Q4','r-',1.2);
%plot_mesh(node,elementV(inactiveElems,:),'Q4','c*');
%plot(gps(:,1),gps(:,2),'+');

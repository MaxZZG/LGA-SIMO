% Data file for a horizontal material interface problem

addpath ../fem_util
addpath ../nurbs-geopdes/inst/
addpath ../nurbs-util/
addpath ../meshing/
addpath ../fem-functions/
addpath ../post-processing/
addpath ../xiga/

%clear all

L = 1;
controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [0;0];
controlPts(1:2,2,1) = [L;0];
controlPts(1:2,1,2) = [0;L];
controlPts(1:2,2,2) = [L;L];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

% h-refinement

refineCount = 3;

for i=1:refineCount
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    newKnotsX = [];
    newKnotsY = [];
    
    % new knots along two directions (uniform)
    
    %newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    
    for it=1:length(uKnotVectorU)-1
        xi1 = uKnotVectorU(it);
        xi2 = uKnotVectorU(it+1);
        le  = (xi2-xi1)/3;
        
        newKnotsX = [newKnotsX xi1+le xi1+2*le];
    end
    
    for it=1:length(uKnotVectorV)-1
        xi1 = uKnotVectorV(it);
        xi2 = uKnotVectorV(it+1);
        le  = (xi2-xi1)/3;
        
        newKnotsY = [newKnotsY xi1+le xi1+2*le];
    end
    
    newKnots  = {newKnotsX newKnotsY};
    
    % h-refinement
    
    solid     = nrbkntins(solid,newKnots);
    
    uKnot      = cell2mat(solid.knots(1));
    vKnot      = cell2mat(solid.knots(2));
end

% and evaluate order 

solid = nrbdegelev(solid,[1 1]); 

% newKnotsX = [0.25 0.55];
% newKnotsY = [0.45 0.55];
% newKnots  = {newKnotsX newKnotsY};
% solid     = nrbkntins(solid,newKnots);

convert2DNurbs

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

%% generate element connectivity ...

generateIGA2DMesh

% plot the mesh

buildVisualizationMesh;

% Geometry of the material interface

linInterfaceGeo   = [0 L/2; L L/2];

% Compute level sets

nNode = size(node,1);
nCPts = size(controlPts,1);
chi   = zeros(1,nNode);
CHI   = zeros(1,nCPts);

x0   = linInterfaceGeo(1,1); y0 = linInterfaceGeo(1,2);
x1   = linInterfaceGeo(2,1); y1 = linInterfaceGeo(2,2);

for i = 1 : nNode
    x                  = node(i,1);
    y                  = node(i,2);
    l                  = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
    phi                = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
    chi(1,i) = phi/l;
end

for i = 1 : nCPts
    x                  = controlPts(i,1);
    y                  = controlPts(i,2);
    l                  = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
    phi                = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
    CHI(1,i) = phi/l;
end

enrich_node = zeros(noCtrPts,1);

count1 = 0;

% loop over elements

for iel = 1 : size(elementV,1)
    sctr    = elementV(iel,:);
    sctrIGA = element(iel,:);
    
    phi  = chi(1,sctr);
        
    if ( max(phi)*min(phi) < 0 )
        count1                 = count1 + 1 ; % ah, one split element
        splitElems(count1)     = iel;                                
        enrich_node(sctrIGA)   = 3;        
    end
end

inc_nodes = find(enrich_node == 3);

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(node,elementV,'Q4','b-',1.2);
% plot the interface
plot(linInterfaceGeo(:,1),linInterfaceGeo(:,2),'k-','Linewidth',1.9);
% plot elements cut by the circle
plot_mesh(node,elementV(splitElems,:),'Q4','r-',1.2);
%plot(gps(:,1),gps(:,2),'+');
n1 = plot(controlPts(inc_nodes,1),controlPts(inc_nodes,2),'r*');
set(n1,'MarkerSize',16,'LineWidth',1.07);

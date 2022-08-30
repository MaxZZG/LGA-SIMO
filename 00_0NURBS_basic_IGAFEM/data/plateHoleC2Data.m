% the following is for the plate with hole
% two element mesh from Hughes et all, CMAME 2005
% k-refinement using the NURBS package
%  (1) order elevation to cubic
%  (2) knot insertion (h-refinement)
% Vinh Phu Nguyen
% Delft University of Technology, The Netherlands

addpath nurbs-geopdes/inst/

noPtsX = 4;
noPtsY = 3;

controlPts = zeros(4,noPtsX,noPtsY);

controlPts(1:2,1,1) = [-1;0];
controlPts(1:2,2,1) = [-1;0.4142135623730951];
controlPts(1:2,3,1) = [-0.4142135623730951; 1];
controlPts(1:2,4,1) = [0;1];

controlPts(1:2,1,2) = [-2.5;0];
controlPts(1:2,2,2) = [-2.5; 0.75];
controlPts(1:2,3,2) = [-0.75; 2.5];
controlPts(1:2,4,2) = [0; 2.5];

controlPts(1:2,1,3) = [-4;0];
controlPts(1:2,2,3) = [-4;4];
controlPts(1:2,3,3) = [-4;4];
controlPts(1:2,4,3) = [0;4];

cont = 0.5*(1+1/sqrt(2));

controlPts(4,:,:)   = 1;
controlPts(4,2,1)   = cont;
controlPts(4,3,1)   = cont;

uKnot = [0 0 0 0.5 1 1 1];
vKnot = [0 0 0 1 1 1];

% build NURBS object

srf = nrbmak(controlPts,{uKnot vKnot});

% order elevate

srf = nrbdegelev(srf,[1 1]); 

% nrbplot(srf,[20 20])
% view([0 90])
% grid off

%%%%%%%%
%% convert NURBS data back to our data structure for analysis

p          = srf.order(1)-1;
q          = srf.order(2)-1;
uKnot      = cell2mat(srf.knots(1));
vKnot      = cell2mat(srf.knots(2));
noPtsX     = length(uKnot)-p-1;
noPtsY     = length(vKnot)-q-1;
weights    = reshape(srf.coefs(4,:,:),noPtsX*noPtsY,1);

controlPts = [];

for ir=1:noPtsY
    controlPts = [controlPts; srf.coefs(1:2,:,ir)'];
end

% refineCount = 1;
% 
% if (refineCount)
%     hRefinement2d
% end
% 
% plotMesh (controlPts,weights,uKnot,vKnot,p,q,50,'r-','try.eps');





       
       

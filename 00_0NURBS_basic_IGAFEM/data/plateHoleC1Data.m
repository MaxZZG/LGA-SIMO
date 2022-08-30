% the following is for the plate with hole
% two element mesh from Hughes et all, CMAME 2005
% Vinh Phu Nguyen
% Delft University of Technology, The Netherlands

controlPts        = zeros(4,3,2);

controlPts(1,1,:) = [-1 0];
controlPts(1,2,:) = [-2.5 0];
controlPts(1,3,:) = [-4 0];

controlPts(2,1,:) = [-1 0.4142135623730951];
controlPts(2,2,:) = [-2.5 0.75];
controlPts(2,3,:) = [-4 4];

controlPts(3,1,:) = [-0.4142135623730951 1];
controlPts(3,2,:) = [-0.75 2.5];
controlPts(3,3,:) = [-4 4];

controlPts(4,1,:) = [0 1];
controlPts(4,2,:) = [0 2.5];
controlPts(4,3,:) = [0 4];

uKnot = [0 0 0 0.5 1 1 1];
vKnot = [0 0 0 1 1 1];

noPtsX = 4;
noPtsY = 3;

controlPts = reshape(controlPts,12,2);

p     = 2;
q     = 2;

cont = 0.5*(1+1/sqrt(2));

weights = [1; cont; cont; 1;
           1; 1; 1; 1;
           1; 1; 1; 1];
       
       

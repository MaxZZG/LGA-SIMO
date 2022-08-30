%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
% 
% L-shaped specimen problem
% 
% Vinh Phu Nguyen, 
% Johns Hopkins University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/

close all
clear all

global p q

refineCount   = 6; % 0: no refinement. Refine mesh with 1, 2 and so on 
noGPs         = 3; % # of Gauss points along one direction


E0  = 1;  % Young's modulus
nu0 = 0.3;  % Poisson's ratio
stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
F = 1;

% COMPUTE ELASTICITY MATRIX

C = elasticityMatrix(E0,nu0,stressState);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LShapedData         %C0 elements
LShapedC1Data        %C1 elements 

% h-refinement here

if (refineCount) 
    hRefinement2d 
end

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

% find boundary nodes for boundary conditions

fixedXNodes  = find(controlPts(:,1)==2*a)';
fixedYNodes  = find(controlPts(:,2)==2*a)';
leftNodes    = find(controlPts(:,1)==0)';

% plot the mesh

plotMesh(controlPts,weights,uKnot,vKnot,p,q,50,'r--','try.eps');

% build connectivity ...

generateIGA2DMesh

leftPoints    = controlPts(leftNodes,:);
leftEdgeMesh  = zeros(noElemsV,q+1);

for i=1:noElemsV
   leftEdgeMesh(i,:) = leftNodes(i:i+q);
end

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix 
f = zeros(noDofs,1);        % external force vector

% essential boundary conditions

uFixed     = zeros(size(fixedXNodes));
vFixed     = zeros(size(fixedYNodes));

udofs      = fixedXNodes;          
vdofs      = fixedYNodes+noCtrPts;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
[W,Q]=quadrature(  noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over elements (knot spans)

for e=1:noElems
   idu    = index(e,1);
   idv    = index(e,2);
   xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
   etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
   
   sctr   = element(e,:);          %  element scatter vector
   sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
   nn     = length(sctr); 
   pts    = controlPts(sctr,:);
   
   B      = zeros(3,2*nn);
 
   % loop over Gauss points 
   
    for gp=1:size(W,1)                        
      pt      = Q(gp,:);                          
      wt      = W(gp);                            
      
      % compute coords in parameter space
      Xi      = parent2ParametricSpace(xiE,pt(1)); 
      Eta     = parent2ParametricSpace(etaE,pt(2)); 
      J2      = jacobianPaPaMapping(xiE,etaE);

      [dRdxi, dRdeta] = NURBS2Dders([Xi; Eta],p,q,uKnot,vKnot,weights');

      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
                  
      jacob = pts'*[dRdxi' dRdeta']; 
      J1    = det(jacob);
      
      % Jacobian inverse and spatial derivatives
            
      dRdx       = [dRdxi' dRdeta']/jacob;
      
      % B matrix
      %        _                                      _
      %        |  N_1,x  N_2,x  ...      0      0  ... |
      %  B  =  |      0      0  ... N_1,y  N_2,y  ... |
      %        |  N_1,y  N_2,y  ... N_1,x  N_2,x  ... |
      %        -                                      -
      
      B(1,1:nn)       = dRdx(:,1)';
      B(2,nn+1:2*nn)  = dRdx(:,2)';
      B(3,1:nn)       = dRdx(:,2)';
      B(3,nn+1:2*nn)  = dRdx(:,1)';
      
      % compute elementary stiffness matrix and
      % assemble it to the global matrix
      
      K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J1 * J2 * wt;
    end
end

% Computing external force

[W1,Q1] = quadrature(2, 'GAUSS', 1 ); 

% Loop over elements along left edge = noElemsV

for e=1:noElemsV
   xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
   conn  = elConnV(e,:);   
   pts   = leftPoints(conn,:);
   sctrx = leftEdgeMesh(e,:);
   
   % loop over Gauss points 
    for gp=1:size(W1,1)                        
      xi      = Q1(gp,:);                          
      wt      = W1(gp);                            
      Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
      J2      = 0.5 * ( xiE(2) - xiE(1) ); 
      
      [N, dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights);
      
      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
      
      x        = N     * pts; % global coord of GP
      jacob1   = dNdxi * pts;
      J1       = norm (jacob1);
         
      f(sctrx) = f(sctrx) + N' * F * J1 * J2 * wt;  
    end
end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
                    % used to keep the  conditioning of the K matrix

applyBC


% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux    = U(1:noCtrPts);
Uy    = U(noCtrPts+1:noDofs);

scale = 0.8;
deformedControlPts = controlPts + scale * [Ux Uy];


uXNode1 = U(1);

% plot stress

vtuFile = '../results/lShaped';

plotStress1










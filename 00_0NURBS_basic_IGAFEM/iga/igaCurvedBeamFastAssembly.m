%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
% 
% Cylinder subject to inner pressure. Only a quarter is modeled.
% 
% Fast assembly using the triplet sparse format. 
%
% Vinh Phu Nguyen, 
% Cardiff University, UK
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../nurbs-geopdes/inst/

close all
clear variables

global p q

noGPs         = 4; % # of Gauss points along one direction


E0          = 3e7;  % Young's modulus
nu0         = 0.25;  % Poisson's ratio
stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS"
pressure    = 3e4; % inner pressure 

% COMPUTE ELASTICITY MATRIX

C = elasticityMatrix(E0,nu0,stressState);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

annularDataGeopdes

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

% find boundary nodes for boundary conditions

fixedXNodes  =  find(controlPts(:,1)==0);
fixedYNodes  =  find(controlPts(:,2)==0);
forcedNodes  =  1:noPtsX:noCtrPts;

% build connectivity ...

generateIGA2DMesh

% build boundary mesh for force vector computation

bndPoints      = controlPts(forcedNodes,:);
rightEdgeMesh  = zeros(noElemsV,q+1);

for i=1:noElemsV
    rightEdgeMesh(i,:) = forcedNodes(i:i+q);
end

% initialization
f = zeros(noDofs,1);        % external force vector

%% assembly
nElNod = size(element,2);
nElDof = nElNod*2;
nElmLK = nElDof^2;
nSprGK = nElmLK*noElems;

jSprLK = 1:nElmLK;
vIdGrd = ones(1,nElDof);

% values, row indices, columns indices of the global K matrix
vSprGK = zeros(nSprGK,1);
jSprRw = zeros(nSprGK,1);
jSprCl = zeros(nSprGK,1);


%% essential boundary conditions

uFixed     = zeros(size(fixedXNodes))';
vFixed     = zeros(size(fixedYNodes))';

udofs      = fixedXNodes;          
vdofs      = fixedYNodes+noCtrPts;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
[W,Q]=quadrature(  noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature

% Assembling system of equations
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
   
   B      = zeros(3,2*nn);
 
   pts    = controlPts(sctr,:);     
   
   Ke = zeros(nElDof,nElDof); % element Ke
   
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
      
      Ke = Ke + B' * C * B * J1 * J2 * wt;
    end
    
    %----------------------------------------------------------------------
    % Global Assembly (Sparse)
    %----------------------------------------------------------------------
    
    vElDof = sctrB';
    mRwGrd = vElDof(:,vIdGrd);
    mClGrd = mRwGrd';
    
    jSprRw(jSprLK) = mRwGrd(:);
    jSprCl(jSprLK) = mClGrd(:);
    vSprGK(jSprLK) = Ke(:);
    
    jSprLK         = jSprLK + nElmLK; % move to the next position 
end

% Here comes the total stiffness matrix in one shot!!!

K = sparse(jSprRw,jSprCl,vSprGK,noDofs,noDofs);

%% Computing external force

[W1,Q1] = quadrature(2, 'GAUSS', 1 ); 

% Loop over elements along left edge = noElemsV

for e=1:noElemsV
   xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
   conn  = elConnV(e,:);
   noFns = length(conn);
   
   sctrx = rightEdgeMesh(e,:);
   sctry = sctrx + noCtrPts;
   
   % loop over Gauss points 
    for gp=1:size(W1,1)                        
      xi      = Q1(gp,:);                          
      wt      = W1(gp);                            
      Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
      J2      = 0.5 * ( xiE(2) - xiE(1) ); 
      
      N       = [];
      dNdxi   = [];
        
      % compute derivative of basis functions w.r.t parameter coord
      
      for in=1:noFns
       [Ni,dNi]  = NURBSbasis (conn(in),q,Xi,vKnot,weights);
       N         = [N Ni];
       dNdxi     = [dNdxi dNi];
      end
      
      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
      
      jacob1   = dNdxi * bndPoints(conn,:);
      J1       = norm (jacob1);
      
      x        = N *bndPoints(conn,:); % global coord of GP
      r        = norm(x);
      Fx       = pressure * x(1,1)/r;
      Fy       = pressure * x(1,2)/r;
      
      f(sctrx) = f(sctrx) + N' * Fx * J1 * J2 * wt;     
      f(sctry) = f(sctry) + N' * Fy * J1 * J2 * wt;  
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

vtuFile     = '../results/cylinderPressure';
plotStress1








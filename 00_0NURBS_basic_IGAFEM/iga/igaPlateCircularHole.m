%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
%
% Infinite plate with a circular hole
%
% Vinh Phu Nguyen,
% Johns Hopkins University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/
addpath ../nurbs-util/

close all
clear all

global p q

refineCount = 3; % 0: no refinement. Refine mesh with 1, 2 and so on

E           = 1e5;  % Young's modulus
nu          = 0.3;  % Poisson’s ratio
stressState = 'PLANE_STRESS';
a           = 1; % hole radius
L           = 4; % length of plate

vtuFile     = '../results/infinitePlateHole';

% Elasticity matrix

C = elasticityMatrix(E,nu,stressState);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plateHoleC1Data
%plateHoleC2Data  % obtained from C1 mesh with k-refinement

% h-refinement here

if (refineCount)
    hRefinement2d
end

noGPs          = p+1; % # of GPs, assume p>=q

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

% Boundary nodes
% 1. Displacement conditions:
%    - On bottom edge (y=0) with u_y = 0
%    - On right  edge (x=0) with u_x = 0

% 2. Traction conditions
%    - On left edge (x=-L) with t = (-sigma_x,-sigma_xy)
%    - On top edge  (y=L) with t = (sigma_xy,sigma_y)

% find boundary nodes for boundary conditions

bottomNodes = find(controlPts(:,2)==0)';
rightNodes  = find(controlPts(:,1)==0)';
leftNodes   = find(controlPts(:,1)==-L)';
topNodes    = find(controlPts(:,2)==L)';

% essential boundary conditions

uFixed     = zeros(size(rightNodes));
vFixed     = zeros(size(bottomNodes));

% plot the mesh

plotMesh(controlPts,weights,uKnot,vKnot,p,q,50,'r-','try.eps');

% build connectivity ...

generateIGA2DMesh

% build a 1D meshes for the right and top edges over which
% a traction is applied

leftPoints    = controlPts(leftNodes,:);
topPoints     = controlPts(topNodes,:);
leftEdgeMesh  = zeros(noElemsV,q+1);
topEdgeMesh   = zeros(noElemsV,q+1);

for i=1:noElemsV
    leftEdgeMesh(i,:) = leftNodes(i:i+q);
end

for i=1:noElemsV    
    topEdgeMesh(i,:) = topNodes(i:i+q);
end

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

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
%%
% Loop over elements (knot spans)

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);         %  element scatter vector
    sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
    nn     = length(sctr);
    B      = zeros(3,2*nn);
    pts    = controlPts(sctr,:);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % compute derivative of basis functions w.r.t parameter coord
        
        [dRdxi, dRdeta] = NURBS2Dders([Xi; Eta],p,q,uKnot,vKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        % Jacobian matrix
        
        jacob     = pts'*[dRdxi' dRdeta'];
        J1        = det(jacob);
        %invJacob  = inv(jacob);
        dRdx      = [dRdxi' dRdeta']/jacob;
        
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
%%
[W1,Q1] = quadrature(8, 'GAUSS', 1 );

% Loop over elements along left edge = noElemsV

for e=1:noElemsV
    xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
    conn  = elConnV(e,:);
    pts   = leftPoints(conn,:);
    sctrx = leftEdgeMesh(e,:);
    sctry = sctrx + noCtrPts;
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights);
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        x        = N     * pts; % global coord of GP
        jacob1   = dNdxi * pts;
        J1       = norm (jacob1);
        
        % compute the exact stresses
        
        str = exact_plate_hole(x,a);
        tx  = -str(1);
        ty  = -str(3);
        fac = J1 * J2 * wt;
        f(sctrx) = f(sctrx) + N' * tx * fac;
        f(sctry) = f(sctry) + N' * ty * fac;
    end
end

% loop over top edge

for e=1:noElemsV
    xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
    conn  = elConnV(e,:);
    pts   = topPoints(conn,:);
    sctrx = topEdgeMesh(e,:);
    sctry = sctrx + noCtrPts;
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N, dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights);
        
        x        = N     * pts; % global coord of GP
        jacob1   = dNdxi * pts;
        J11      = norm (jacob1);
        
        % compute the exact stresses
        
        str = exact_plate_hole(x,a);
        tx  = str(3) ;
        ty  = str(2);
        fac = J11 * J2 * wt;
        f(sctrx) = f(sctrx) + N' * tx * fac;
        f(sctry) = f(sctry) + N' * ty * fac;
    end
end

% penalty method to enforce the two overlapping control points
% at the top left corner to have the same displacements.
% u_a = u_b: see IFEM lecture note, C. Fellipa, Colorado.

w     = 1000;
sctr  = [leftNodes(end-1) leftNodes(end)];
sctry = sctr + noCtrPts;
penaltyStiffness = w*[1 -1;-1 1];
K(sctr,sctr)   = K(sctr,sctr)   + penaltyStiffness;
K(sctry,sctry) = K(sctry,sctry) + penaltyStiffness;

%%
disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

udofs = rightNodes;            % global indecies  of the fixed x disps
vdofs = bottomNodes+noCtrPts;  % global indecies  of the fixed y disps

applyBC

% SOLVE SYSTEM

disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux    = U(1:noCtrPts);
Uy    = U(noCtrPts+1:noDofs);

%%%%%%%%%%%%%%%%%%%%%
%plot stress field

plotStress1

% compute energy and displacement norms

computeNorms

% compare exact and numerical displacements
ux_exact = zeros(1, size(node,1));
uy_exact = zeros(1, size(node,1));

for i = 1 : size(node,1)
    x           = node(i,:) ;
    [exact_stress, exact_disp] = exact_solution_hole(x,a,E,nu);
    ux_exact(i) = exact_disp(1);
    uy_exact(i) = exact_disp(2);
end
% ----------------------------------

% --------------------------------------------
% Plot both exact and numerical deformed shape
fac=600;
figure
hold on
h = plot(node(:,1)+fac*dispX,node(:,2)+fac*dispY,'rs');
set(h,'MarkerSize',7);
h = plot(node(:,1)+fac*ux_exact',node(:,2)+fac*uy_exact','b*');
set(h,'MarkerSize',7);
title('Exact and numerical deformed shape')
legend('XIGA','Exact')
axis equal





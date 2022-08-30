%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended Isogeometric analysis for two dimensional linear elastic
% fracture mechanics problems.
%
% Center crack infinite plate in mode II loading
% Dirichlet BCs = exact displacements imposed with penalty method.
%
% Vinh Phu Nguyen,
% Johns Hopkins University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/
addpath ../nurbs-geopdes/inst/
addpath ../nurbs-util/

close all
clear all

global p q controlPts weights element xCrack xTips crack_node
global uKnot vKnot noDofs levelSets split_nodes spit_elem noGPs1
global jDomainFac
jDomainFac = 6;

E0          = 1e7;  % Young modulus
nu0         = 0.3;  % Poisson ratio
stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
sigmato     = 1e4;

% COMPUTE ELASTICITY MATRIX
if ( strcmp(stressState,'PLANE_STRESS')  )      % Plane Strain case
    C=E0/(1-nu0^2)*[  1      nu0          0;
        nu0        1          0;
        0        0  (1-nu0)/2  ];
else                                            % Plane Strain case
    C=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0      nu0        0;
        nu0    1-nu0        0;
        0        0  1/2-nu0 ];
end

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edgeInfiniteCenterCrackC0Data
%edgeInfiniteCenterCrackC2Data

noGPs         = p+1; % # of Gauss points along one direction
noGPs1        = p+1;

% find boundary nodes for boundary conditions

bottomNodes =  find(controlPts(:,2)==0); 
rightNodes  =  find(controlPts(:,1)==D); 
topNodes    =  find(controlPts(:,2)==D); 
leftNodes   =  find(controlPts(:,1)==0); 

% build boundary mesh 

bottomPoints    = controlPts(bottomNodes,:);
rightPoints     = controlPts(rightNodes,:);
topPoints       = controlPts(topNodes,:);
leftPoints      = controlPts(leftNodes,:);

bottomEdgeMesh  = zeros(noElemsU,p+1);
rightEdgeMesh   = zeros(noElemsV,q+1);
topEdgeMesh     = zeros(noElemsU,p+1);
leftEdgeMesh    = zeros(noElemsU,q+1);

for i=1:noElemsU
    bottomEdgeMesh(i,:) = bottomNodes(i:i+p);
    topEdgeMesh(i,:)    = topNodes(i:i+p);
end

for i=1:noElemsV
    rightEdgeMesh(i,:) = rightNodes(i:i+q);    
    leftEdgeMesh(i,:)  = leftNodes(i:i+q);  
end

% Each split node is enriched by ONE function, H(x)
% Each tip node is enriched by FOUR functions, B_i(x)
% then, the total dofs is :
% total dof = numnode*nsd + numsplitnode*1*nsd + numstipnode*4*nsd
% here, two dimension, nsd = 2

noDofs = noDofs + size(split_nodes,1)*1*2 + ...
    size(tip_nodes,1)*4*2;

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

% We use fictitious nodes/control points to handle these
% additional dofs. At a H(x) enriched node, we add one fantom node and
% at tip enriched node, four fantom nodes are added. These fictitious nodes
% are numbered from the total number of true nodes, ie, from numnode+1 ...

pos    = zeros(noCtrPts,1);
nsnode = 0 ;
ntnode = 0 ;

for i = 1 : noCtrPts
    if (enrich_node(i) == 1)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4) + 1 ;
        nsnode = nsnode + 1 ;
    elseif (enrich_node(i) == 2)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4) + 1 ;
        ntnode = ntnode + 1 ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    nn     = length(sctr);
    
    % -----------------------------------------------
    % Choose Gauss quadrature rules for elements
    
    if     (ismember(e,split_elem))     % split element
        [W,Q] = quadrature(12,'GAUSS',2);
    elseif (ismember(e,tip_elem))       % tip element
        [W,Q] = quadrature(12,'GAUSS',2);
    elseif ( any(intersect(tip_nodes,sctr)) ~= 0)% having tip enriched nodes
        [W,Q] = quadrature(12,'GAUSS',2);
    else
        [W,Q] = quadrature(noGPs,'GAUSS',2);
    end
    
    % Determine the position in the global matrix K
    
    sctrB = assembly(e,enrich_node,pos);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
      
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
    
        % compute derivative of basis functions w.r.t parameter coord

        [N, dRdxi, dRdeta] = NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');

        % B matrix
        
        [B,J1] = BMatrixXIGA(Xi,Eta,e,enrich_node,N,dRdxi,dRdeta);
        
        % Stiffness matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B'*C*B*W(gp)*J1*J2;   
    end
end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

% If penalty method is use, then one must modify the stiffness matrix and
% the nodal force vector
% $K_{ij}$ = Kij - alpha \int phi_i phi_j d \gamma_u
% fj       = fj  - alpha \int phi_i u_bar d \gamma_u

% Loop over elements along top edge 

fu = zeros(noDofs,1);
k  = zeros(noDofs,noDofs);

[ff, kk] = penaltyBoundaryCondition(elRangeU,elConnU,noElemsU,...
    bottomEdgeMesh,bottomPoints,E0,nu0,stressState,...
    sigmato,xTip,seg,cracklength,pos,2);

fu   = fu + ff;
k    = k + kk;

[ff, kk] = penaltyBoundaryCondition(elRangeU,elConnU,noElemsU,...
    topEdgeMesh,topPoints,E0,nu0,stressState,...
    sigmato,xTip,seg,cracklength,pos,2);

fu   = fu + ff;
k    = k + kk;

[ff, kk, en] = penaltyBoundaryCondition(elRangeV,elConnV,noElemsV,...
    rightEdgeMesh,rightPoints,E0,nu0,stressState,...
    sigmato,xTip,seg,cracklength,pos,2);

fu   = fu + ff;
k    = k + kk;

[W1,Q1] = quadrature(8, 'GAUSS', 1 );

% Loop over elements along left edge = noElemsV
% left edge: natural BCs

for e=1:noElemsV
    xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
    conn  = elConnV(e,:);
    noFns = length(conn);
    sctrx = 2*leftEdgeMesh(e,:)-1;
    sctry = 2*leftEdgeMesh(e,:);
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );

        [N, dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights);
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        x        = N    *leftPoints(conn,:); % global coord of GP
        jacob1   = dNdxi*leftPoints(conn,:);
        J1       = norm (jacob1);
        
        [sigmaxx,sigmayy,sigmaxy] = exactStressModeII(x,E0,nu0,...
           sigmato,xTip,seg,cracklength);
        
        tx = -sigmaxx;
        ty = -sigmaxy;
       
        f(sctrx) = f(sctrx) + N' * tx * J1 * J2 * wt;
        f(sctry) = f(sctry) + N' * ty * J1 * J2 * wt;
    end
end


% modified system of equations

%alpha = 0.05*max(max(C)) ;           % penalty number
alpha = 1e10;
f = f - alpha*fu;
K = K - alpha*k;


% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux = U(1:2:noDofs);
Uy = U(2:2:noDofs);

% compute stresses and displacements at mesh vertices
disp([num2str(toc),'  COMPUTING STRESSES '])

vtuFile      = '../results/infiniteCrackModeII';
vtuCrackFile = '../results/infiniteCrackModeII-cracked';

plotStressXIGA;

% plot stress
fac=1e6;

figure
clf
hold on
plot_field(node+fac*[dispX dispY],elementV,'Q4',sigmaXY);
hold on
colorbar
title('Sigma_{xy}')
axis off
plot_mesh(node+fac*[dispX dispY],elementV,'Q4','g.-',1.2);

figure
clf
plot_field(node+fac*[dispX dispY],elementV,'Q4',displacement(:,:,1));
plot_mesh(node+fac*[dispX dispY],elementV,'Q4','g.-',1.2);
hold on
colorbar
title('Displacement in x direction')
axis off
set(gcf,'color','white')

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

computeInteractionIntegrals

% Compute SIFs from I integral

Knum = I.*E0/(2*(1-nu0^2)) % plain strain 

% Compute the exact SIFs
 
Kexact = sigmato*sqrt(pi*cracklength)


%%%%%%%%%%%%
eSigmaXX = zeros(1, size(node,1));
eSigmaYY = zeros(1, size(node,1));
eSigmaXY = zeros(1, size(node,1));
eSigmaVM = zeros(1, size(node,1));

ux_exact = zeros(1, size(node,1));
uy_exact = zeros(1, size(node,1));

for i = 1 : size(node,1)
    x = node(i,:) ;
    [ux1,uy1] = exactDispModeII(x,E0,nu0,stressState,sigmato,xTip,seg,cracklength);
    ux_exact(i) = ux1;
    uy_exact(i) = uy1;
    [sigmaxx,sigmayy,sigmaxy] = exactStressModeII(x,E0,nu0,...
                                  sigmato,xTip,seg,cracklength);
    eSigmaXX(i) = sigmaxx;
    eSigmaYY(i) = sigmayy;
    eSigmaXY(i) = sigmaxy;
    % von Mises stress
    eSigmaVM(i) = sqrt(sigmaxx^2+sigmayy^2-sigmaxx*sigmayy+3*sigmaxy^2);
end

vtuFile = '../results/infiniteCrackModeIIExact';
VTKPostProcess(node,elementV,'Quad4',vtuFile,...
    [eSigmaXX' eSigmaYY' eSigmaXY' eSigmaVM'],[ux_exact' uy_exact'])

% --------------------------------------------
% Plot both exact and numerical deformed shape

figure
hold on
h = plot(node(:,1)+fac*dispX,node(:,2)+fac*dispY,'rs');
set(h,'MarkerSize',7);
h = plot(node(:,1)+fac*ux_exact',node(:,2)+fac*uy_exact','b*');
set(h,'MarkerSize',7);
title('Exact and numerical deformed shape')
legend('XIGA','Exact')
axis equal
axis off
set(gcf,'color','white')
% --------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post-processing for cracks
crackedMeshNURBS


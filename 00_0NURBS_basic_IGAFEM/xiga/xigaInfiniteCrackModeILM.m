%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended Isogeometric analysis for two dimensional linear elastic
% fracture mechanics problems.
%
% Center crack infinite plate in tension
% Dirichlet BCs = exact displacements imposed with the
% Lagrange multiplier method.
%
% Vinh Phu Nguyen,
% Delft University of Technology
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
global uKnot vKnot noDofs levelSets split_nodes spit_elem
global jDomainFac

jDomainFac  = 2;

E0          = 1e7;  % Young modulus
nu0         = 0.3;  % Poisson ratio
stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
sigmato     = 1e4;

vtuFile      = '../results/infiniteCrackModeILM';
vtuCrackFile = '../results/infiniteCrackModeILM-cracked';

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
%edgeInfiniteCenterCrackC4Data
%edgeInfiniteCenterCrackC6Data

noGPs         = p+1; % # of Gauss points along one direction
noGP1         = p+1;

% find boundary nodes for boundary conditions

bottomNodes =  find(node(:,2)==0);
rightNodes  =  find(node(:,1)==D);
topNodes    =  find(node(:,2)==D);
leftNodes   =  find(node(:,1)==0);

bottomNodesIGA =  find(controlPts(:,2)==0);
rightNodesIGA  =  find(controlPts(:,1)==D);
topNodesIGA    =  find(controlPts(:,2)==D);
leftNodesIGA   =  find(controlPts(:,1)==0);

disp_nodes = [bottomNodesIGA ; rightNodesIGA; topNodesIGA];
disp_nodes = unique(disp_nodes);
num_disp_nodes = length(disp_nodes);

topNodesIGA    = sort(topNodesIGA,'descend');

% build boundary mesh

bottomPoints    = controlPts(bottomNodesIGA,:);
rightPoints     = controlPts(rightNodesIGA,:);
topPoints       = controlPts(topNodesIGA,:);
leftPoints      = controlPts(leftNodesIGA,:);

% linear two node line elements 

bottomEdgeMesh  = zeros(length(bottomNodesIGA)-1,2);
rightEdgeMesh   = zeros(length(rightNodesIGA)-1,2);
topEdgeMesh     = zeros(length(topNodesIGA)-1,2);
leftEdgeMesh    = zeros(length(leftNodesIGA)-1,2);

bottomEdgeMeshIGA  = zeros(noElemsU,p+1);
rightEdgeMeshIGA   = zeros(noElemsV,q+1);
topEdgeMeshIGA     = zeros(noElemsU,p+1);
leftEdgeMeshIGA    = zeros(noElemsU,q+1);

for i=1:noElemsU
    bottomEdgeMeshIGA(i,:) = bottomNodesIGA(i:i+p);    
    topEdgeMeshIGA(i,:)    = topNodesIGA(i:i+p);
end

for i=1:noElemsV    
    rightEdgeMeshIGA(i,:) = rightNodesIGA(i:i+p);    
    leftEdgeMeshIGA(i,:)  = leftNodesIGA(i:i+p);
end

for i=1:length(bottomNodesIGA)-1
    bottomEdgeMesh(i,:) = bottomNodesIGA(i:i+1);    
    topEdgeMesh(i,:)    = topNodesIGA(i:i+1);    
end

for i=1:length(rightNodesIGA)-1
    rightEdgeMesh(i,:) = rightNodesIGA(i:i+1);    
    leftEdgeMesh(i,:)  = leftNodesIGA(i:i+1);    
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
    elseif (ismember(e,tip_elem))   % tip element
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

[W1,Q1] = quadrature(5, 'GAUSS', 1 );

qk = zeros(1,2*num_disp_nodes);
G  = zeros(noDofs,2*num_disp_nodes);
    
% Loop over L2 finite elements not NURBS elements along 
% bottom edge, right edge and top edge in THAT ORDER !!!
% to ensure the correct assembly of Lagrange multiplier unknowns.
% Due to this, the top edge is reserved see the SORT command 
% in line 74.

m1 = 0 ;

for i = 1 : length(bottomEdgeMesh)
    sctr = bottomEdgeMesh(i,:);
    m1   = m1 + 1 ;
    m2   = m1 + 1 ;
    pts  = controlPts(sctr,:);  
    
    for gp = 1:size(W1,1)
        pt        = Q1(gp,:);
        wt        = W1(gp);
        [N,dNdxi] = lagrange_basis('L2',pt);       
        J0        = dNdxi'*pts;
        detJ      = norm(J0) ;
        pt        = N' * pts; % global GP
        
        % compute exact displacement
        
        [ux,uy] = exactDispModeI(pt,E0,nu0,stressState,sigmato,xTip,seg,...
                                    cracklength);
        
        N1 = N(1) ; N2 = N(2) ;
        
        % qk vector
        fac1       = wt * detJ * N1;
        fac2       = wt * detJ * N2;

        qk(2*m1-1) = qk(2*m1-1) - fac1 * ux;
        qk(2*m1)   = qk(2*m1)   - fac1 * uy;
        qk(2*m2-1) = qk(2*m2-1) - fac2 * ux;
        qk(2*m2)   = qk(2*m2)   - fac2 * uy;
        
        % G matrix
        
        xi    = inverseMapping1DNURBS (pt(1),uKnot,bottomPoints,p,1,elConnU,1);
        uspan = FindSpan(noPtsX-1,p,xi,uKnot);
        nurbsCon = bottomNodesIGA(uspan-p+1:uspan+1);
        [N, dNdxi] = NURBS1DBasisDers(xi,p,uKnot,weights);
                
        for j = 1 : length(nurbsCon)
            row1 = 2*nurbsCon(j)-1 ;
            row2 = 2*nurbsCon(j)   ;
            G1   = - wt * detJ * N(j) * [N1 0 ; 0 N1];
            G2   = - wt * detJ * N(j) * [N2 0 ; 0 N2];
            G(row1:row2,2*m1-1:2*m1) = G(row1:row2,2*m1-1:2*m1) + G1;
            G(row1:row2,2*m2-1:2*m2) = G(row1:row2,2*m2-1:2*m2) + G2;
        end
    end
end


for i = 1 : length(rightEdgeMesh)
    sctr = rightEdgeMesh(i,:);
    m1   = m1 + 1 ;
    m2   = m1 + 1 ;
    pts  = controlPts(sctr,:);  

    for gp = 1:size(W1,1)
        pt        = Q1(gp,:);
        wt        = W1(gp);
        [N,dNdxi] = lagrange_basis('L2',pt);       
        J0        = dNdxi'*pts;
        detJ      = norm(J0) ;
        pt        = N'*pts; % global GP
        
        % compute exact displacement
        
        [ux,uy] = exactDispModeI(pt,E0,nu0,stressState,sigmato,xTip,seg,...
                                    cracklength);
        
        N1 = N(1) ; N2 = N(2) ;
        
        % qk vector
        fac1 = wt * detJ * N1;
        fac2 = wt * detJ * N2;

        qk(2*m1-1) = qk(2*m1-1) - fac1 * ux;
        qk(2*m1)   = qk(2*m1)   - fac1 * uy;
        qk(2*m2-1) = qk(2*m2-1) - fac2 * ux;
        qk(2*m2)   = qk(2*m2)   - fac2 * uy;
        
        % G matrix
        
        xi        = inverseMapping1DNURBS (pt(2),vKnot,rightPoints,q,2,elConnV,1);
        uspan     = FindSpan(noPtsY-1,q,xi,vKnot);
        nurbsCon  = rightNodesIGA(uspan-q+1:uspan+1);
        [N, dNdxi] = NURBS1DBasisDers(xi,q,vKnot,weights);
                
        for j = 1 : length(nurbsCon)
            row1 = 2*nurbsCon(j)-1 ;
            row2 = 2*nurbsCon(j)   ;
            G1   = - wt * detJ * N(j) * [N1 0 ; 0 N1];
            G2   = - wt * detJ * N(j) * [N2 0 ; 0 N2];
            G(row1:row2,2*m1-1:2*m1) = G(row1:row2,2*m1-1:2*m1) + G1;
            G(row1:row2,2*m2-1:2*m2) = G(row1:row2,2*m2-1:2*m2) + G2;
        end
    end
end

for i = 1 : length(topEdgeMesh)
    sctr = topEdgeMesh(i,:);
    m1   = m1 + 1 ;
    m2   = m1 + 1 ;    
    pts  = controlPts(sctr,:);
        
    for gp = 1:size(W1,1)
        pt        = Q1(gp,:);
        wt        = W1(gp);
        [N,dNdxi] = lagrange_basis('L2',pt);       
        J0        = dNdxi'*pts;
        detJ      = norm(J0) ;
        pt        = N' * pts; % global GP
        
        % compute exact displacement
        
        [ux,uy] = exactDispModeI(pt,E0,nu0,stressState,sigmato,xTip,seg,...
                                    cracklength);
        
        N1 = N(1) ; N2 = N(2) ;
        
        % qk vector
        fac1       = wt * detJ * N1;
        fac2       = wt * detJ * N2;

        qk(2*m1-1) = qk(2*m1-1) - fac1 * ux;
        qk(2*m1)   = qk(2*m1)   - fac1 * uy;
        qk(2*m2-1) = qk(2*m2-1) - fac2 * ux;
        qk(2*m2)   = qk(2*m2)   - fac2 * uy;
        
        % G matrix
        
        xi        = inverseMapping1DNURBS (pt(1),uKnot,topPoints,p,1,elConnU,-1);
        uspan     = FindSpan(noPtsX-1,p,xi,uKnot);
        nurbsCon  = topNodesIGA(uspan-p+1:uspan+1);
        [N, dNdxi] = NURBS1DBasisDers(xi,p,uKnot,weights);
                
        for j = 1 : length(nurbsCon)
            row1 = 2*nurbsCon(j)-1 ;
            row2 = 2*nurbsCon(j)   ;
            G1   = - wt * detJ * N(j) * [N1 0 ; 0 N1];
            G2   = - wt * detJ * N(j) * [N2 0 ; 0 N2];
            G(row1:row2,2*m1-1:2*m1) = G(row1:row2,2*m1-1:2*m1) + G1;
            G(row1:row2,2*m2-1:2*m2) = G(row1:row2,2*m2-1:2*m2) + G2;
        end
    end
end

% Loop over elements along left edge = noElemsV
% left edge: natural BCs

for e=1:noElemsV
    xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
    conn  = elConnV(e,:);
    sctrx = 2*leftEdgeMeshIGA(e,:)-1;
    sctry = 2*leftEdgeMeshIGA(e,:);
    pts   = leftPoints(conn,:);

    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N, dNdxi] = NURBS1DBasisDers(Xi,p,uKnot,weights);
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        x        = N    *pts; % global coord of GP
        jacob1   = dNdxi*pts;
        J1       = norm (jacob1);
        
	    % compute exact stresses

        [sigmaxx,sigmayy,sigmaxy] = exactStressModeI(x,E0,nu0,...
                                   sigmato,xTip,seg,cracklength);
        
	    % then the tractions, note that normal vector is [-1,0]    

        tx = -sigmaxx;
        ty = -sigmaxy;
        
        f(sctrx) = f(sctrx) + N' * tx * J1 * J2 * wt;
        f(sctry) = f(sctry) + N' * ty * J1 * J2 * wt;
    end
end


% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])

f = [f;qk'];                                % f = {f;qk}
m = ([K G; G' zeros(num_disp_nodes*2)]);    % m = [K GG;GG' 0]
d = m\f;                                    % d = {u;lamda}

% just get nodal parameters u_i, not need Lagrange multipliers
U = d(1:noDofs);
clear d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux = U(1:2:noDofs);
Uy = U(2:2:noDofs);

% compute stresses and displacements at mesh vertices
disp([num2str(toc),'  COMPUTING STRESSES '])

plotStressXIGA;

% plot stress

fac=1e6;
stressComp=2;
figure
clf
hold on
plot_field(node+fac*[dispX dispY],elementV,'Q4',sigmaYY);
hold on
colorbar
title('Stress in x direction')
axis off
plot_mesh(node+fac*[dispX dispY],elementV,'Q4','g.-',1.2);
set(gcf,'color','white') 

figure
clf
plot_field(node+fac*[dispX dispY],elementV,'Q4',displacement(:,:,2));
plot_mesh(node+fac*[dispX dispY],elementV,'Q4','g.-',1.2);
hold on
colorbar
title('Displacement in y direction')
axis off
set(gcf,'color','white')   

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([num2str(toc),'      SIFs COMPUTATION'])

computeInteractionIntegrals

% Compute SIFs from I integral

Knum = I.*E0/(2*(1-nu0^2)) % plain strain

% Compute the exact SIFs

Kexact = sigmato*sqrt(pi*cracklength)


%%%%%%%%%%%%
ux_exact = zeros(1, size(node,1));
uy_exact = zeros(1, size(node,1));

for i = 1 : size(node,1)
    x = node(i,:) ;
    [ux1,uy1] = exactDispModeI(x,E0,nu0,stressState,sigmato,xTip,seg,cracklength) ;
    ux_exact(i) = ux1;
    uy_exact(i) = uy1;
end
% ----------------------------------

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
% --------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post-processing for cracks

crackedMeshNURBS


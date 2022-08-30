%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended Isogeometric analysis for two dimensional linear elastic
% fracture mechanics problems.
%
% Center crack infinite plate in tension
% Dirichlet BCs = exact displacements imposed with penalty method.
%
% Vinh Phu Nguyen,
% Delft University of Technology, The Netherlands
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path(path,'~/Downloads/codeIGA_crack/fem_util');
path(path,' ../C_files/');
path(path,' ../data/');
path(path,' ../meshing/');
path(path,' ../post-processing/');
path(path,' ../fem-functions/');
path(path,' ../analytical-solutions/');
path(path,' ../nurbs-geopdes/inst/');
path(path,' ../nurbs-util/');
path(path,' ../integration/');
path(path,' ./function');

close all
clear all

global p q controlPts weights element xCrack xTips crack_node
global uKnot vKnot noDofs levelSets split_nodes spit_elem noGPs1
global cont node

E0          = 1e7;  % Young modulus
nu0         = 0.3;  % Poisson ratio
stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
sigmato     = 1e4;

vtuFile      = 'infinieCrackModeI';
vtuCrackFile = 'infinieCrackModeI-cracked';

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

%edgeInfiniteCenterCrackC0Data
edgeCrack2D_data

noGPs         = p+1; % # of Gauss points along one direction
noGPs1        = p+1;

% find boundary nodes for boundary conditions

% bottomNodes =  find(controlPts(:,2)==0);
% rightNodes  =  find(controlPts(:,1)==L);
% topNodes    =  find(controlPts(:,2)==D);
% leftNodes   =  find(controlPts(:,1)==0);

bottomNodes = 1 : noPtsX;
topNodes    = noPtsX*(noPtsY-1)+1 : noCtrPts;
topNodes    = sort(topNodes,'descend');

rightNodes  = noPtsX : noPtsX : noCtrPts;
leftNodes   = 1: noPtsX : noCtrPts;
% build boundary mesh

bottomPoints    = controlPts(bottomNodes,:);
rightPoints     = controlPts(rightNodes,:);
topPoints       = controlPts(topNodes,:);
leftPoints      = controlPts(leftNodes,:);

% Loop over elements along top edge 

[elRangeU,elConnU] = buildConnectivity(p,uKnot,Numx);
bottomEdgeMesh  = bottomNodes(elConnU);
topEdgeMesh  = topNodes(elConnU);

[elRangeV,elConnV] = buildConnectivity(q,vKnot,Numy);
rightEdgeMesh  = rightNodes(elConnV);
leftEdgeMesh  = leftNodes(elConnV);


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

%%% Assembling system of equation
%%% Stiffness matrix and external force vector

disp([num2str(toc,'%0.6g'),'  ASSEMBLING THE SYSTEM'])

% Loop over elements (knot spans)
tol=1e-8;
%
for e=1:numelem
    ni = inn(ien(e,1),1);% get NURBS coordinates
    nj = inn(ien(e,1),2);
    if (abs(uKnot(ni)-uKnot(ni+1))>tol)&&(abs(vKnot(nj)-vKnot(nj+1))>tol)
        xiE=[uKnot(ni), uKnot(ni+1)];
        etaE=[vKnot(nj), vKnot(nj+1)];
               
        sctr   = element(e,:);          %  element scatter vector
        sctrV  = elementV(e,:);         %  element scatter vector Q4
        pts    = controlPts(sctr,:);
        nn     = length(sctr);
   
    % -----------------------------------------------
    % Choose Gauss quadrature rules for elements
%         [W,Q] = gaussPointForEnrichedElement(e,noGPs,sctrV,levelSets,xTip,split_elem,...
%                                         tip_elem,tip_nodes);
        
    
    if     (ismember(e,split_elem)) % split element, 13 GPs/subtriangle
%         [W,Q] = discontQ4quad(7,levelSets(1,elementV(e,:),1));
        [W,Q] = quadrature(20,'GAUSS',2);
    elseif (ismember(e,tip_elem))   % tip element
%         [W,Q] = disTipQ4quad(7,levelSets(1,elementV(e,:)),node(sctrV,:),xTip);
        [W,Q] = quadrature(20,'GAUSS',2);
    elseif ( any(intersect(tip_nodes,sctr)) ~= 0)% having tip enriched nodes
        [W,Q] = quadrature(12,'GAUSS',2);
    else
        [W,Q] = quadrature(10,'GAUSS',2);
    end
    
    % Determine the position in the global matrix K
    
        sctrB = assembly(e,element,enrich_node,pos);
    
    % loop over Gauss points
    
        for gp=1:size(W,1)
            pt      = Q(gp,:);
        % compute coords in parameter space
            Xi   = ((xiE(2) - xiE(1))*pt(1)+xiE(2) + xiE(1))/2;
            Eta  = ((etaE(2)-etaE(1))*pt(2)+etaE(2)+etaE(1))/2;
            J2   = (xiE(2) - xiE(1))*(etaE(2)-etaE(1))/4;
            
%         Xi      = parent2ParametricSpace(xiE,pt(1));
%         Eta     = parent2ParametricSpace(etaE,pt(2));
%         J2      = jacobianPaPaMapping(xiE,etaE);
    
        % compute derivative of basis functions w.r.t parameter coord
        
            [N dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');
         
        % B matrix 
        
            [B,J1] = BMatrixXIGA(Xi,Eta,e,enrich_node,N,dRdxi,dRdeta);
            
        
        % Stiffness matrix
             
            K(sctrB,sctrB) = K(sctrB,sctrB) + B'*C*B*W(gp)*J1*J2;
        end

    end
end
%}
noElemsU=Numx;
noElemsV=Numy;
noElems =numelem;

%% Boundary conditions

disp([num2str(toc,4),'  APPLYING DIRICHLET BOUNDARY CONDITIONS'])

% If penalty method is use, then one must modify the stiffness matrix and
% the nodal force vector
% $K_{ij}$ = Kij - alpha \int phi_i phi_j d \gamma_u
% fj       = fj  - alpha \int phi_i u_bar d \gamma_u

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Least square method for determining the control
% points displacements

noBndElems = Numx*2 + Numy;
bndElement = zeros(noBndElems,p+1);
dispNodes   = [bottomNodes , rightNodes(2:end), ...
               topNodes(2:end)];
noDispNodes = length(dispNodes);

for ie=1:Numx
    bndElement(ie,:) = [ie:ie+p]+(ie-1)*(Rep-1);
end

start = max(max(bndElement));

for ie=1:Numy
    bndElement(ie+Numx,:) = [start:start+p]+(ie-1)*(Rep-1);
    start = start + 1;
end

start = max(max(bndElement));

for ie=1:Numx
    bndElement(ie+Numx+Numy,:) = [start:start+p]+(ie-1)*(Rep-1);
    start = start + 1;
end

A  = zeros(noDispNodes,noDispNodes);
bx = zeros(noDispNodes,1);
by = zeros(noDispNodes,1);

% bottom edge

noxC   = 9;

for ie=1:Numx
    sctr   = bottomEdgeMesh(ie,:);
    pts    = controlPts(sctr,:);
    sctrA  = bndElement(ie,:);
    xiE    = elRangeU(ie,:);
    xiArr  = linspace(xiE(1),xiE(2),noxC);
   
    for ic=2:noxC-1                
        xi = xiArr(ic);
        [N dNdxi] = NURBS1DBasisDers(xi,p,uKnot,weights);
        
        A(sctrA,sctrA) = A(sctrA,sctrA) + N'*N;
        
        x        = N    *pts; 
        
        % exact displacements
        
        [ux,uy] = exactDispModeI(x,E0,nu0,stressState,sigmato,xTip,seg,...
            cracklength);
        
        bx(sctrA) = bx(sctrA) + ux*N';
        by(sctrA) = by(sctrA) + uy*N';
    end
end

% right edge

for ie=1:Numy
    sctr   = rightEdgeMesh(ie,:);
    pts    = controlPts(sctr,:);
    sctrA  = bndElement(ie+Numx,:);
    xiE    = elRangeV(ie,:);
    xiArr  = linspace(xiE(1),xiE(2),noxC);
    
    for ic=2:noxC-1                     
        xi             = xiArr(ic);        
        [N dNdxi]      = NURBS1DBasisDers(xi,q,vKnot,weights);        
        A(sctrA,sctrA) = A(sctrA,sctrA) + N'*N;        
        x              = N*pts; 
        
        % exact displacements
        
        [ux,uy] = exactDispModeI(x,E0,nu0,stressState,...
                                 sigmato,xTip,seg,cracklength);
        
        bx(sctrA) = bx(sctrA) + ux*N';
        by(sctrA) = by(sctrA) + uy*N';
    end
end

% top edge

for ie=1:Numx
    sctr   = topEdgeMesh(ie,:);
    pts    = controlPts(sctr,:);    
    sctrA  = bndElement(ie+Numx+Numy,:);
    xiE    = elRangeU(ie,:);
    xiArr  = linspace(xiE(1),xiE(2),noxC);
        
    for ic=2:noxC-1                        
        xi        = xiArr(ic);        
        [N dNdxi] = NURBS1DBasisDers(xi,p,uKnot,weights);        
        A(sctrA,sctrA) = A(sctrA,sctrA) + N'*N;        
        x         = N*pts; 
        
        % exact displacements
        
        [ux,uy] = exactDispModeI(x,E0,nu0,stressState,...
                                 sigmato,xTip,seg,cracklength);
        
        bx(sctrA) = bx(sctrA) + ux*N';
        by(sctrA) = by(sctrA) + uy*N';
    end
end

% solve the system Aq_x=bx and Aq_y=by
qx=A\bx; 
qy=A\by;

% [LL UU] = lu(A);
% qxTemp  = LL\bx;
% qyTemp  = LL\by;
% qx      = UU\qxTemp;
% qy      = UU\qyTemp;

%%%%

uFixed     = qx';
vFixed     = qy';

udofs      = 2*dispNodes-1;
vdofs      = 2*dispNodes;


clear A bx by


disp([num2str(toc),'  APPLYING NEUMANN BOUNDARY  CONDITIONS'])

[W1,Q1] = quadrature(noGPs1+1, 'GAUSS', 1 );

% Loop over elements along left edge = noElemsV
% left edge: natural BCs

for e=1:Numy
    xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
%     conn  = elConnV(e,:);
    sctrx = 2*leftEdgeMesh(e,:)-1;
    sctry = 2*leftEdgeMesh(e,:);
%     pts   = leftPoints(conn,:);
    pts   = controlPts(leftEdgeMesh(e,:),:);
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N dNdxi] = NURBS1DBasisDers(Xi,p,uKnot,weights);
        
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


%{
% Loop over elements along top edge

fu = zeros(noDofs,1);
k  = zeros(noDofs,noDofs);

[ff kk en] = penaltyBoundaryCondition(elRangeU,elConnU,Numx,...
    bottomEdgeMesh,bottomPoints,E0,nu0,stressState,...
    sigmato,xTip,seg,cracklength,pos,1);

fu   = fu + ff;
k    = k + kk;

[ff kk en] = penaltyBoundaryCondition(elRangeU,elConnU,Numx,...
    topEdgeMesh,topPoints,E0,nu0,stressState,...
    sigmato,xTip,seg,cracklength,pos,1);

fu   = fu + ff;
k    = k + kk;

[ff kk en] = penaltyBoundaryCondition(elRangeV,elConnV,Numy,...
    rightEdgeMesh,rightPoints,E0,nu0,stressState,...
    sigmato,xTip,seg,cracklength,pos,1);

fu   = fu + ff;
k    = k + kk;
% 
% clear kk ff
[W1,Q1] = quadrature(noGPs1, 'GAUSS', 1 );

%% left edge: natural BCs
% Loop over elements along left edge = noElemsV

for e=1:Numy
    xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
    conn  = elConnV(e,:);
    sctrx = 2*leftEdgeMesh(e,:)-1;
    sctry = 2*leftEdgeMesh(e,:);
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N dNdxi] = NURBS1DBasisDers(Xi,p,uKnot,weights);
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        x        = N    *leftPoints(conn,:); % global coord of GP
        jacob1   = dNdxi*leftPoints(conn,:);
        J1       = norm (jacob1);
        
        [sigmaxx,sigmayy,sigmaxy] = exactStressModeI(x,E0,nu0,...
            sigmato,xTip,seg,cracklength);
        
        tx = -sigmaxx;
        ty = -sigmaxy;
        
        f(sctrx) = f(sctrx) + N' * tx * J1 * J2 * wt;
        f(sctry) = f(sctry) + N' * ty * J1 * J2 * wt;
    end
end


% take the enrichment into account when compute the
% external force vector

% for e=1:noElemsV
%     sctr = leftEdgeMesh(e,:);
%     le   = length(sctr);
%
%     if (length(intersect(sctr,split_nodes))~=le)
%         continue
%     end
%
%     xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
%     conn  = elConnV(e,:);
%     noFns = length(conn);
%     sctrx = 2*pos(sctr)-1;
%     sctry = 2*pos(sctr);
%
%     xCr    = reshape(xCrack(1,:,:),2,2);
%     Ne     = zeros(le,1);
%
%     % loop over Gauss points
%     for gp=1:size(W1,1)
%         xi      = Q1(gp,:);
%         wt      = W1(gp);
%         Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
%         J2      = 0.5 * ( xiE(2) - xiE(1) );
%
%         N       = [];
%         dNdxi   = [];
%
%         % compute derivative of basis functions w.r.t parameter coord
%
%         for in=1:noFns
%             [Ni,dNi]  = NURBSbasis (conn(in),q,Xi,vKnot,weights);
%             N         = [N Ni];
%             dNdxi     = [dNdxi dNi];
%         end
%
%         % compute the jacobian of physical and parameter domain mapping
%         % then the derivative w.r.t spatial physical coordinates
%
%         x        = N    *leftPoints(conn,:); % global coord of GP
%         jacob1   = dNdxi*leftPoints(conn,:);
%         J1       = norm (jacob1);
%
%         [sigmaxx,sigmayy,sigmaxy] = exactStressModeI(x,E0,nu0,...
%                                    sigmato,xTip,seg,cracklength);
%
%         tx = -sigmaxx;
%         ty = -sigmaxy;
%
%         % Enrichment function, H(x) at global Gauss point
%         dist = signed_distance(xCr,x);
%         Hgp  = heaviside(dist);
%
%         % Enrichment function, H(x) at node "in"
%
%         for i=1:length(sctr)
%             nodeId = sctr(i);
%             dist   = signed_distance(xCr,controlPts(nodeId,:));
%             H      = heaviside(dist);H=0;
%             Ne(i)  = N(i)*(Hgp-H);
%         end
%
%         f(sctrx) = f(sctrx) + Ne * tx * J1 * J2 * wt;
%         f(sctry) = f(sctry) + Ne * ty * J1 * J2 * wt;
%     end
% end

% modified system of equations

% alpha = 0.05*max(max(C)) ;           % penalty number
% alpha = 1e10;
% f = f - alpha*fu;
% K = K - alpha*k;

% clear fu k
%}
% SOLVE SYSTEM
disp([num2str(toc,'%0.6g'),'  SOLVING THE SYSTEM'])

applyBC

U=K\f;
% break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ux = U(1:2:noDofs);
% Uy = U(2:2:noDofs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
disp([num2str(toc,'%0.6g'),'  SIFs COMPUTATION'])

computeInteractionIntegrals

% Compute SIFs from I integral

% Knum = I.*E0/(2*(1-nu0^2)) % plain strain

% Compute the exact SIFs

Kexact = sigmato*sqrt(pi*cracklength)

% compute stresses and displacements at mesh vertices
disp([num2str(toc,'%0.6g'), '  COMPUTING STRESSES '])

% % plot stress through Gauss Points
% plotStressGP   

% plot stress through element nodes
plotStress;

%
%%%%%%%%%%%%

for i = 1 : size(node,1)
    x = node(i,:) ;
    [ux1,uy1] = exactDispModeI(x,E0,nu0,stressState,sigmato,xTip,seg,cracklength);
    ux_exact(i) = ux1;
    uy_exact(i) = uy1;
    [sigmaxx,sigmayy,sigmaxy] = exactStressModeI(x,E0,nu0,...
        sigmato,xTip,seg,cracklength);
    eSigmaXX(i) = sigmaxx;
    eSigmaYY(i) = sigmayy;
    eSigmaXY(i) = sigmaxy;
    % von Mises stress
    eSigmaVM(i) = sqrt(sigmaxx^2+sigmayy^2-sigmaxx*sigmayy+3*sigmaxy^2);
end

vtuFile = 'infiniteCrackModeIExact';
VTKPostProcess(node,elementV,2,'Quad4',vtuFile,...
    [eSigmaXX' eSigmaYY' eSigmaXY' eSigmaVM'],[ux_exact' uy_exact'])

% ----------------------------------

% --------------------------------------------
% Plot both exact and numerical deformed shape
fac=3000000;
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

% crackedMeshNURBS

%}

%{
% % ======================= ERROR NORM ================================
% calculate error norm
disp([num2str(toc,'%0.6g'),'  ERROR NORM COMPUTATION'])
dnorm=0; enorm=0;

for e=1:numelem
    ni = inn(ien(e,1),1);% get NURBS coordinates
    nj = inn(ien(e,1),2);
    if (abs(uKnot(ni)-uKnot(ni+1))>tol)&&(abs(vKnot(nj)-vKnot(nj+1))>tol)
        xiE=[uKnot(ni), uKnot(ni+1)];
        etaE=[vKnot(nj), vKnot(nj+1)];
        
        sctr   = element(e,:);          %  element scatter vector
        nn     = length(sctr);
        pts    = controlPts(sctr,:);
    
        uspan = FindSpan(noPtsX-1,p,xiE(1),uKnot);
        vspan = FindSpan(noPtsY-1,q,etaE(1),vKnot);
    
        elemDisp  = element_disp(e,pos,enrich_node,U);
        
        
    if     (ismember(e,split_elem)) % split element, 13 GPs/subtriangle
%         [W,Q] = discontQ4quad(7,levelSets(1,elementV(e,:),1));
        [W,Q] = quadrature(20,'GAUSS',2);
    elseif (ismember(e,tip_elem))   % tip element
%         [W,Q] = disTipQ4quad(7,levelSets(1,elementV(e,:)),node(sctrV,:),xTip);
        [W,Q] = quadrature(20,'GAUSS',2);
    elseif ( any(intersect(tip_nodes,sctr)) ~= 0)% having tip enriched nodes
        [W,Q] = quadrature(10,'GAUSS',2);
    else
        [W,Q] = quadrature(noGPs1,'GAUSS',2);
    end
    
    
    % loop over Gauss points
         for gp=1:size(W,1)
            pt      = Q(gp,:);
        % compute coords in parameter space
            Xi   = ((xiE(2) - xiE(1))*pt(1)+xiE(2) + xiE(1))/2;
            Eta  = ((etaE(2)-etaE(1))*pt(2)+etaE(2)+etaE(1))/2;
            J2   = (xiE(2) - xiE(1))*(etaE(2)-etaE(1))/4;
            [N dRdxi dRdeta] = NURBS2DBasisDersSpecial([Xi; Eta],...
                               p,q,uKnot,vKnot,weights',[uspan;vspan]);
            
            [B,J1] = BMatrixXIGA(Xi,Eta,e,enrich_node,N, dRdxi, dRdeta);
            [exN]  = NMatrixXIGA(e,enrich_node,N);
            
            strain  = B*elemDisp;
            stress  = C*strain;
            
            uh      = exN*[elemDisp(1:2:end)  elemDisp(2:2:end)];
            
            x0      = N*controlPts(sctr,:);
            
            [ux1,uy1] = exactDispModeI(x0,E0,nu0,stressState,sigmato,xTip,seg,cracklength);
            [sigmaxx,sigmayy,sigmaxy] = exactStressModeI(x0,E0,nu0,sigmato,xTip,seg,cracklength);
            
            dnorm=dnorm+([ux1,uy1]-uh)*([ux1,uy1]-uh)'*W(gp)*J1*J2;
            enorm=enorm+(stress'-[sigmaxx,sigmayy,sigmaxy])*inv(C)* ...
                  (stress'-[sigmaxx,sigmayy,sigmaxy])'*W(gp)*J1*J2;            
         end
    end
end

dnorm=sqrt(dnorm)  % error norm displacement
enorm=sqrt(enorm)  % error norm energy
    
%}    
    
    
    

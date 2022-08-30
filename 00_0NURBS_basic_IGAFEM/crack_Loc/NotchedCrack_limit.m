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

path(path,'../fem_util');
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


% path (path, './functions_nurbs')
%path (path, '~/mosek/6/tools/platform/linux32x86/bin')
%path (path, '~/mosek/6/toolbox/r2009b')


close all
clear all

global p q controlPts weights element xCrack xTips crack_node
global uKnot vKnot noDofs levelSets split_nodes spit_elem noGPs1
global cont node

E0          = 1e4;  % Young modulus
nu0         = 0.3;  % Poisson ratio
stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
sigmap      = 1;    % Yield stress MPa

% Loading
sigma = 1  ;

vtuFile      = 'infinieCrackModeI';
vtuCrackFile = 'infinieCrackModeI-cracked';

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%edgeInfiniteCenterCrackC0Data
NotchedCrack_limit_data

noGPs         = p; % # of Gauss points along one direction
noGPs1        = p;

% find boundary nodes for boundary conditions

% bottomNodes =  find(controlPts(:,2)==0);
% rightNodes  =  find(controlPts(:,1)==L);
% topNodes    =  find(controlPts(:,2)==D);
% leftNodes   =  find(controlPts(:,1)==0);

bottomNodes = 1 : noPtsX;
topNodes    = noPtsX*(noPtsY-1)+1 : noCtrPts;
% topNodes    = sort(topNodes,'descend');

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
% u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector
S = 0;                      % count the number of GaussPoints
coordGP = [];               % gcoordinate of Gauss Points

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

%          [W,Q] = quadrature(noGPs,'GAUSS',2);
         
        if     (ismember(e,split_elem)) % split element, 13 GPs/subtriangle
%                 [W,Q] = discontQ4quad(7,levelSets(1,elementV(e,:),1));
               [W,Q] = quadrature(2*p,'GAUSS',2);
        elseif (ismember(e,tip_elem))   % tip element
%                  [W,Q] = disTipQ4quad(7,levelSets(1,elementV(e,:)),node(sctrV,:),xTip);
              [W,Q] = quadrature(2*p,'GAUSS',2);
        elseif ( any(intersect(tip_nodes,sctr)) ~= 0)% having tip enriched nodes
                [W,Q] = quadrature(p+1,'GAUSS',2);
        else
                [W,Q] = quadrature(noGPs,'GAUSS',2);
        end
    
    % Determine the position in the global matrix K
    
        sctrB = assembly(e,element,enrich_node,pos);
    
    % loop over Gauss points
    
        for gp=1:size(W,1)
            pt   = Q(gp,:);
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

            S  =  S + 1;
            
            Gpnt = N * pts;               % global GP
            coordGP = [coordGP;Gpnt];
            
            BB = zeros (3,noDofs);
            BB(:, sctrB) = B;
            Bstrain{S}   = BB;
            WeJ(S)       = W(gp)*J1*J2;
            clear BB
        
%         %Stiffness matrix             
%             K(sctrB,sctrB) = K(sctrB,sctrB) + B'*C*B*W(gp)*J1*J2;
        end     % end loop for GPs

    end  % end if
end % end loop for elements

% ************************
%    NODAL FORCE VECTOR
% ************************
disp([num2str(toc),'   NODAL FORCE VECTOR COMPUTATION'])

[W1,Q1] = quadrature(noGPs1, 'GAUSS', 1 );

% The right edge is applied a traction along Y direction
for e=1:Numy   
    xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
    sctry = 2*rightEdgeMesh(e,:)-1;
    pts   = controlPts(rightEdgeMesh(e,:),:);
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
       [N dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights);        
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob1   = dNdxi * pts;
        J1       = norm (jacob1);
              
        f(sctry) = f(sctry) + N' * (sigma) * J1 * J2 * wt;
    end
end

%{
% % BC
fixedXNodes= leftNodes';
fixedYNodes= 1;
uFixed     = zeros(size(fixedXNodes))';
vFixed     = zeros(size(fixedYNodes))';

udofs      = 2*fixedXNodes-1;
vdofs      = 2*fixedYNodes;


applyBC

% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noElems=numelem;
plotStress;

fac=300;
figure
hold on
h = plot(node(:,1)+fac*dispX,node(:,2)+fac*dispY,'rs');
set(h,'MarkerSize',7);
% h = plot(node(:,1)+fac*ux_exact',node(:,2)+fac*uy_exact','b*');
% set(h,'MarkerSize',7);
title('Exact and numerical deformed shape')
legend('XIGA','Exact')
axis equal
axis off
set(gcf,'color','white')
break
%}

% The left edge is applied a traction along Y direction
for e=1:Numy   
    xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
    sctry = 2*leftEdgeMesh(e,:)-1;
    pts   = controlPts(leftEdgeMesh(e,:),:);
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
       [N dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights);        
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob1   = dNdxi * pts;
        J1       = norm (jacob1);
              
        f(sctry) = f(sctry) + N' * (-sigma) * J1 * J2 * wt;
    end
end

Wex = f'; % External Work
clear f

% % % boundary condition
% fixedEle = Numx * (Numy+1)/2;
% fixedXNodes = find(controlPts(element(fixedEle,:),1)==1); %bottom edge
% fixedXNodes = element(fixedEle, fixedXNodes);
% bcdof = 2.*fixedXNodes - 1;

% -------------------------------------

% **********************************
%    ESSENTIAL BOUNDARY CONDITION
% **********************************
disp([num2str(toc),'   Build matrices for OP'])

% +----------------------------+
% | OBJECTIVE FUNCTION AMPHA+  |
% +----------------------------+
k1 = S; % nel*nG for t_i __ --- sodinhtai * socanh = cac gia tri ti
%k1 = ned;
k2 = 3*k1; % 3*k1 for r_i: r1,r2,r

%sdof = total_unknown; %cac bien chuyen vi mo rong cua tat ca cac dinh tai
var = noDofs + k1 + k2;
f = zeros(var,1); % [u1 v1...un vn ...bien mo rong t11...t(nG+nel)r1...rk2 mp1....mp_nel]
%(tong bac tu do + tong diem Gauss + 3*tong diem Gauss)
for i = 1:k1
    f(noDofs+i,1) = sigmap*WeJ(i);
end
% External work of unitary and boundary conditions
% Wex = force'; % External Work
% Aeq = [Wex];

% % boundary condition
% sbc = size(bcdof,2);
% Bc = zeros(sbc,noDofs);
% for i = 1:sbc
%     loca = bcdof(i);
%     Bc(i,loca) = 1;
% end
% Aeq = [Wex; Bc];
% beq = zeros(1+sbc,1);

Aeq = [Wex];
beq(1,1) = 1; % external work of unitary

a1 = size(Aeq,1);
% trans to X
A1 = sparse([[Aeq] [zeros(a1,k1+k2)]]);
% these variables appear when assigning C'k = r_i
% [A2] = added_variables_EXFEM(noDofs,k1,k2,Bstrain,S);
[A2] = added_variables_strain_ex(noDofs,k1,k2,Bstrain,S);

aeq = [[A1]; [A2]];
b = [[beq];[zeros(k2,1)]];

% cones
for i = 1:k1
    co(i,:) = [noDofs + i, noDofs+k1+3*i-2, noDofs+k1+3*i-1, noDofs+k1+3*i];
end

% +--------------+
% | OPTIMIZATION |
% +--------------+
% Mosek optimisation tool: MOSEKOPT
disp([num2str(toc),'   OPTIMIZATION PROCESS'])

'solving'
clear prob
clear param
prob.c = f';
prob.a = aeq;
prob.buc = b';
prob.blc = b';
%lx = [[-inf*ones(nno,1)]; [zeros(nno,1)]; [-inf*ones(3*nno,1)]];
prob.blx = [];
prob.bux = [];
% define cones
prob.cell = cell(k1,1);
for i = 1:k1
    prob.cones{i}.type = 'MSK_CT_QUAD';
    prob.cones{i}.sub  = co(i,:);
    %prob.cones{i}.sub  = [nno + i, 2*nno + 3*i-2, 2*nno + 3*i-1, 2*nno + 3*i]
end

param =[];
param.MSK_IPAR_PRESOLVE_USE = 'MSK_OFF';
[r,res] = mosekopt('minimize',prob,param);
try
    % Display the optimal solution .
    u = res.sol.itr.xx;
    result = f'*u/sigmap;
catch
    fprintf ('MSKERROR : Could not get solution')
end

% ------------------------ plot displacement ------------------------------

var
%Nel
result
% break
figure

scale = 0.15;
dis_u = u(1:2:2*noCtrPts-1);
dis_v = u(2:2:2*noCtrPts);
gcoord2(:,1) = controlPts(:,1)+ dis_u*scale;
gcoord2(:,2) = controlPts(:,2)+ dis_v*scale;
gcoord2(:,3) = weights;
icount=0;
for j = 1:noPtsY
      for i = 1:noPtsX
          icount = icount + 1;
            control(i,j,:) = gcoord2(icount,:);
      end
end  

CP1(:,:,1)=control(:,:,1);
CP1(:,:,2)=control(:,:,2);
CP1(:,:,3)=ones(size(control(:,1,1),1),size(control(1,:,1),2));
CP1(:,:,4)=control(:,:,3);

% nodes = connecElementQ4(noPtsX-1,noPtsY-1);
% patch('faces',nodes,'vertices',controlPts,'facecolor','g','edgecolor','none');
axis equal, hold on
%plotNURBS_surf(u_knot,v_knot,CP); hold on
plotNURBS_surf(uKnot,vKnot,CP1); axis equal 
% axis ([-0.5 5.5 0 5.5]);



diss = u(noDofs+1:noDofs+k1);
% get decompused triangles
tri = delaunay(coordGP(:,1),coordGP(:,2));
tri = tricheck(coordGP,tri);

figure
clf
plot_field(coordGP,tri,'T3',diss);
hold on
colorbar
title('Displacement in y direction')
axis off; view(2); hold on;
plot (xCr1(:,1),xCr1(:,2),'-w', 'LineWidth', 2)    
plot (xCr2(:,1),xCr2(:,2),'-w', 'LineWidth', 2)        

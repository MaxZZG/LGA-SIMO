clear all
deg=2; ref=1;  %%% refinement

%%%%%%%% section hole plate
CP=zeros(7,6,4);
p=2;
q=2;

U = [0 0 0 1/5 2/5 3/5 4/5 1 1 1]; %n=7
V = [0 0 0 1/4 2/4 3/4 1 1 1]; %m=6

R=1;

CP(:,:,1)=[0 0 0 4*R 4*R 4*R;
           0 0.32*R 0.74*R 3.25*R 3.67*R 4*R;
           0.5*R 1.44*R 1.70*R 2.25*R 2.55*R 3.5*R;
           1.25*R 1.50*R 1.75*R 2.25*R 2.50*R 2.75*R;
           0.5*R 1.44*R 1.70*R 2.25*R 2.55*R 3.5*R;
           0 0.32*R 0.74*R 3.25*R 3.67*R 4*R;
           0 0 0 4*R 4*R 4*R];
       
CP(:,:,2)=[R 4*R 4*R 4*R 4*R R;
           R 3.49*R 3.45*R 3.46*R 3.49*R R;
           R 0.77*R 0.74*R 0.75*R 0.78*R R; 
           0 0 0 0 0 0;
           -R -0.77*R -0.74*R -0.75*R -0.78*R -R;
           -R -3.49*R -3.45*R -3.46*R -3.49*R -R;
           -R -4*R -4*R -4*R -4*R -R];
       
CP(:,:,3)=ones(7,6);

CP(:,:,4)=[1 1 1 1 1 1; 1 1 1 1 1 1;
          2/3 1 1 1 1 2/3; 
          4/9 1 1 1 1 4/9;
          2/3 1 1 1 1 2/3;
          1 1 1 1 1 1; 1 1 1 1 1 1];



 %=====================================================================
% REFINE
[CP,U,V,p,q] = degree_elevate_surf_repeated(p,q,U,V,CP,deg-p,deg-q);

plotNURBS_surf(U,V,CP); hold on
plot_ctrlnet(CP,'ro');view(2)
axis equal
break
R1 = refinement_vec_repeated(U,ref,p-1);
R2 = refinement_vec_repeated(V,ref,p-1);

[CP,u_knot,v_knot] = knot_refine_surf(p,q,U,V,CP,R1,R2);

plotNURBS_surf(u_knot,v_knot,CP); hold on
plot_ctrlnet(CP,'ro');view(2)

% section computational IgA
p=deg;
q=p;
ngauss=p;                 % number of gauss point in integration

% input control net
b_net(:,:,1)=CP(:,:,1);
b_net(:,:,2)=CP(:,:,2);
b_net(:,:,3)=CP(:,:,4);

mcp = length(b_net(:,1,1));
ncp = length(b_net(1,:,1));

clear CP U V R1 R2

gcoord(:,1)=reshape(b_net(:,:,1),mcp*ncp,1);
gcoord(:,2)=reshape(b_net(:,:,2),mcp*ncp,1);
gcoord(:,3)=reshape(b_net(:,:,3),mcp*ncp,1);

nodes=connecElementQ4(mcp-1,ncp-1);
%{
for iel=1:size(nodes,1)
    nd=nodes(iel,:);
    x=gcoord(nd,1);
    y=gcoord(nd,2);    
%    xx=mean(x); 
%    yy=mean(y); 
     patch(x,y,'w')
     text(x(1),y(1),num2str(nd(1)));
     text(x(2),y(2),num2str(nd(2)));
     text(x(3),y(3),num2str(nd(3)));
     text(x(4),y(4),num2str(nd(4)));
end

%}

%number of elements, degrees of freedom
nnode=mcp*ncp;         % number of control point
nshl = (p+1)*(q+1);    % number of local shape functions
numx=mcp-p;
numy=ncp-q;
nel=numx*numy; % number of element
nsd=2;                       % number of spatial dimension
ndof=2;
sdof=nnode*ndof;
nnel=nshl;

for i=1:deg+1
   leftEdge(:,i) = [i:1:mcp-p-1+i]';
   rightEdge(:,i) = [mcp*(ncp-1)+i:1:mcp*ncp-p-1+i]';
   upperEdge(:,i) = [1+(i-1)*mcp:mcp:mcp*(ncp-p-1+i)]';
   lowerEdge(:,i) = [(1+i-1)*mcp:mcp:mcp*(ncp-p-1+i)]';
end


%%for force boundary

left=upperEdge(1:ref+p-2,:);
between=upperEdge(ref+p-2+1:size(upperEdge,1)-ref-p+2,:);
right=upperEdge(size(upperEdge,1)-ref-1:size(upperEdge,1),:);


% % for sym boundary condition
lower_sym=[mcp:mcp:mcp*ncp];
%left_sym=[1:mcp:mcp*ncp-mcp+1];
bcdof = [lower_sym*2];


%----------------------------------------------
%  initialization of matrices and vectors
%----------------------------------------------
K=sparse(sdof,sdof); % stiffness matrix
F=zeros(sdof,1);
S = 0; % variable to calculate number of nGauss

% %Read in mesh information from file mesh.dat to load b_net
% icount=0;
% for j = 1:ncp
%       for i = 1:mcp
%           icount = icount + 1;
%             b_net(i,j,:) = gcoord(icount,:);
%       end
% end  

% Read in Master/Slave relationships. These arise from periodic boundary
%conditions, corners in the geometry, or any other case where two
%control points must always remain equal. iper(nnode)
iper=1:mcp*ncp;

%Read in Boundary/Face information: closed_u_flag, closed_v_flag
% flag determines if surface is closed
%  in some direction
closed_u_flag= 0;
closed_v_flag= 0;
%Set number of faces
nedge = 2*(ncp-q)*(1-closed_u_flag)+2*(mcp-p)*(1-closed_v_flag) ;
%Read in boundary condition indicator at each node, in each direction
% 0 - nothing, 1 - Dir., 2 - Neu., 3 - Periodic
%  ibc(nnode,nsd)
ibc=zeros(size(gcoord,1),2);
%ibc(left_sym,1)=1;
ibc(lower_sym,2)=1;

%Read in boundary condition indicator on each face, in each direction
%0 - nothing, 1 - constraint, 2 - load
%IBC_FACE(NFACE,NSD)

ibc_edge=[zeros(size(leftEdge,1),2);zeros(size(lowerEdge,1),2); zeros(size(rightEdge,1),2); zeros(size(left,1),2);[zeros(size(between,1),1) 2*ones(size(between,1),1)];zeros(size(right,1),2)];

% generate connectivities     
[ien,inn]=genIEN_INN_2D(p,q,mcp,ncp);
%generate mesh face/interior loading attributes 
[edge_ien,edge_or]=genEdge_IEN_2D(ien,closed_u_flag,closed_v_flag,nedge);

%*************************************************************************
%Loop over interior elements to build lhsk and force vector contributions 
%of elements (body loads).
%*************************************************************************
[Bstrain,WeJ,S]=KmatNurbs2D_limit(ngauss,nel,inn,ien,b_net);
%[Bstrain,WeJ,S]=KmatNurbs2D_limit_new(ngauss,nel,inn,ien,b_net);

% Loop over boundary edges to get force vector contributions
F=EdgeForces_2D_hole(ngauss,nedge,inn,edge_ien,edge_or,ibc_edge,ibc,b_net,F,load1,load2);

force = F;% External Work
clear F
% -------------------------------------

% **********************************
%    ESSENTIAL BOUNDARY CONDITION
% **********************************

'Build matrices for OP'
% +----------------------------+
% | OBJECTIVE FUNCTION AMPHA+  |
% +----------------------------+
k1 = S; % nel*nG for t_i __ --- sodinhtai * socanh = cac gia tri ti
%k1 = ned;
k2 = 3*k1; % 3*k1 for r_i: r1,r2,r

%sdof = total_unknown; %cac bien chuyen vi mo rong cua tat ca cac dinh tai
var = sdof + k1 + k2;
f = zeros(var,1); % [u1 v1...un vn ...bien mo rong t11...t(nG+nel)r1...rk2 mp1....mp_nel]
%(tong bac tu do + tong diem Gauss + 3*tong diem Gauss)
for i = 1:k1
    f(sdof+i,1) = sigmap*WeJ(i);
end
% External work of unitary and boundary conditions
Wex = force'; % External Work
%Aeq = [Wex];
% boundary condition
sbc = size(bcdof,2);
Bc = zeros(sbc,sdof);
for i = 1:sbc
    loca = bcdof(i);
    Bc(i,loca) = 1;
end
%
Aeq = sparse([Wex; Bc]);
beq = zeros(1+sbc,1);
beq(1,1) = 1; % external work of unitary


a1 = size(Aeq,1);
% trans to X
A1 = sparse([[Aeq] [zeros(a1,k1+k2)]]);
% these variables appear when assigning C'k = r_i
[A2] = added_variables_EXFEM(sdof,k1,k2,Bstrain,S);

aeq = [[A1]; [A2]];
b = [[beq];[zeros(k2,1)]];

% cones
for i = 1:k1
    co(i,:) = [sdof + i, sdof+k1+3*i-2, sdof+k1+3*i-1, sdof+k1+3*i];
end
% +--------------+
% | OPTIMIZATION |
% +--------------+
% Mosek optimisation tool: MOSEKOPT
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

nel
%Nel
result
figure

scale = 0.25;
dis_u = u(1:2:2*nnode-1);
dis_v = u(2:2:2*nnode);
gcoord2(:,1) = gcoord(:,1)+ dis_u*scale;
gcoord2(:,2) = gcoord(:,2)+ dis_v*scale;
gcoord2(:,3) = gcoord(:,3)+ scale;
icount=0;
for j = 1:ncp
      for i = 1:mcp
          icount = icount + 1;
            control(i,j,:) = gcoord2(icount,:);
      end
end  

CP1(:,:,1)=control(:,:,1);
CP1(:,:,2)=control(:,:,2);
CP1(:,:,3)=ones(size(control(:,1,1),1),size(control(1,:,1),2));
CP1(:,:,4)=control(:,:,3);

patch('faces',nodes,'vertices',gcoord,'facecolor','g','edgecolor','none');
axis equal, hold on
%plotNURBS_surf(u_knot,v_knot,CP); hold on
plotNURBS_surf(u_knot,v_knot,CP1); axis equal 
axis ([-0.5 4.5 0 4.5]);

figure
hold on

diss = u(sdof+1:sdof+k1);

%
% plotDissipation(nodes,gcoord,diss); hold on
% break
for i = 1:size(nodes,1)
    X=gcoord(nodes(i,:),1);  Y=gcoord(nodes(i,:),2);
    kCell = convhull(X,Y);
    XX = X(kCell(1:end-1)); YY = Y(kCell(1:end-1));
    fill(XX,YY,diss(i),'edgecolor','none'); 
end
axis equal off


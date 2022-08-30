addpath ../fem_util
addpath ../gmshFiles/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/

clear all
colordef black
state = 0;

% ******************************************************************************
% ***                             I N P U T                                  ***
% ******************************************************************************
tic;


% MATERIAL PROPERTIES
E0  = 30e6;  % Young?s modulus
nu0 = 0.3;  % Poisson?s ratio
% BEAM PROPERTIES
L  = 46;     % length of the beam
c  = 6;      % the distance of the outer fiber of the beam from the mid-line

% MESH PROPERTIES
elemType = 'Q4'; % the element type used in the FEM simulation;
numy     = 4;
numx     = 18;
plotMesh = 1;

% three node constant strain triangular element, ?Q4? is for
% a four node quadrilateral element, and ?Q9? is for a nine
% node quadrilateral element.
% the number of elements in the x-direction (beam length)
% and in the y-direciton.
% A flag that if set to 1 plots the initial mesh (to make sure
% that the mesh is correct)

% TIP LOAD
P = 1000; % the peak magnitude of the traction at the right edge
% STRESS ASSUMPTION


% ******************************************************************************
% ***                    P R E - P R O C E S S I N G                         ***
% ******************************************************************************
I0=2*c^3/3;     % the second polar moment of inertia of the beam cross-section.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE ELASTICITY MATRIX
C=E0/(1-nu0^2)*[   1      nu0          0;
    nu0        1          0;
    0        0  (1-nu0)/2 ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE FINITE ELEMENT MESH
%

% Here we gnerate the finte element mesh (using the approriate elements).
% I won?t go into too much detail about how to use these functions.  If
% one is interested one can type - help ?function name? at the matlab comand
% line to find out more about it.
%
% The folowing data structures are used to describe the finite element
% discretization:
%
%   node    - is a matrix of the node coordinates, i.e. node(I,j) -> x_Ij
%   element - is a matrix of element connectivities, i.e. the connectivity
%             of element e is given by > element(e,:) -> [n1 n2 n3 ...];
%
% To apply boundary conditions a description of the boundaries is needed.  To
% accomplish this we use a separate finite element discretization for each
% boundary.  For a 2D problem the boundary discretization is a set of 1D elements.
%
%   rightEdge -  a element connectivity matrix for the right edge
%   leftEdge  -  I?ll give you three guesses
%
% These connectivity matricies refer to the node numbers defined in the
% coordinate matrix node.
disp([num2str(toc),'   GENERATING MESH'])
switch elemType
    case 'Q4'           % here we generate the mesh of Q4 elements
        nnx=numx+1;
        nny=numy+1;
        node=square_node_array([0 -c],[L -c],[L c],[0 c],nnx,nny);
        inc_u=1;
        inc_v=nnx;
        node_pattern=[ 1 2 nnx+2 nnx+1 ];
        element=make_elem(node_pattern,numx,numy,inc_u,inc_v);
    case 'Q9'           % here we generate a mehs of Q9 elements
        nnx=2*numx+1;
        nny=2*numy+1;
        node=square_node_array([0 -c],[L -c],[L c],[0 c],nnx,nny);
        inc_u=2;
        inc_v=2*nnx;
        node_pattern=[ 1 3 2*nnx+3 2*nnx+1 2 nnx+3 2*nnx+2 nnx+1 nnx+2 ];
        element=make_elem(node_pattern,numx,numy,inc_u,inc_v);
    otherwise %?T3?    % and last but not least T3 elements
        nnx=numx+1;
        nny=numy+1;
        node=square_node_array([0 -c],[L -c],[L c],[0 c],nnx,nny);
        node_pattern1=[ 1 2 nnx+1 ];
        node_pattern2=[ 2 nnx+2 nnx+1 ];
        inc_u=1;
        inc_v=nnx;
        element=[make_elem(node_pattern1,numx,numy,inc_u,inc_v);
            make_elem(node_pattern2,numx,numy,inc_u,inc_v) ];
end
% DEFINE BOUNDARIES

%   Here we define the boundary discretizations.
uln=nnx*(nny-1)+1;
urn=nnx*nny;
lrn=nnx;
lln=1;
cln=nnx*(nny-1)/2+1;  % node number at (0,0)
switch elemType
    
    % upper left node number
    % upper right node number
    % lower right node number
    % lower left node number
    case 'Q9'
        rightEdge=[ lrn:2*nnx:(uln-1); (lrn+2*nnx):2*nnx:urn; (lrn+nnx):2*nnx:urn ]';
        leftEdge =[ uln:-2*nnx:(lrn+1); (uln-2*nnx):-2*nnx:1; (uln-nnx):-2*nnx:1 ]';
        edgeElemType='L3';
    otherwise  % same discretizations for Q4 and T3 meshes
        rightEdge=[ lrn:nnx:(uln-1); (lrn+nnx):nnx:urn ]';
        leftEdge =[ uln:-nnx:(lrn+1); (uln-nnx):-nnx:1 ]';
        edgeElemType='L2';
end
% GET NODES ON DISPLACEMENT BOUNDARY
%      Here we get the nodes on the essential boundaries

fixedNode=unique(leftEdge);  % a vector of the node numbers which are fixed in
% the x direction

uFixed=zeros(1,length(fixedNode))';     % a vector of the x-displacement for the nodes
vFixed=zeros(1,length(fixedNode))';     % and the y-displacements for fixedNodeY


numnode=size(node,1);    % number of nodes
numelem=size(element,1); % number of elements
% PLOT MESH
if ( plotMesh )  % if plotMesh==1 we will plot the mesh
    clf
    plot_mesh(node,element,elemType,'g.-');
    hold on
    plot_mesh(node,rightEdge,edgeElemType,'bo-');
    plot_mesh(node,leftEdge,edgeElemType,'bo-');
    plot(node(fixedNode,1),node(fixedNode,2),'r>');
    axis off
    axis([0 L -c c])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES
%
% Here we define the system data structures
%   f - is the nodal force vector.  It?s structure is the same as U,
%
% K-
%
disp([num2str(toc),' INITIALIZING DATA STRUCTURES'])
U=zeros(2*numnode,1); % nodal displacement vector
f=zeros(2*numnode,1);          % external load vector
K=sparse(2*numnode,2*numnode); % stiffness matrix
% a vector of indicies that quickly address the x and y portions of the data
% strtuctures so U(xs) returns U_x the nodal x-displacements
xs=1:numnode;                  % x portion of u and v vectors
ys=(numnode+1):2*numnode;      % y portion of u and v vectors
% ******************************************************************************
% ***                          P R O C E S S I N G                           ***
% ******************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE EXTERNAL FORCES
%   integrate the tractions on the left and right edges
disp([num2str(toc),'   COMPUTING EXTERNAL LOADS'])
switch elemType  % define quadrature rule
    case 'Q9'
        [W,Q]=quadrature( 4, 'GAUSS', 1 ); %  four point quadrature
    otherwise
        [W,Q]=quadrature( 3, 'GAUSS', 1 ); % three point quadrature
end

% RIGHT EDGE
for e=1:size(rightEdge,1) % loop over the elements in the right edge
    sctr=rightEdge(e,:);  % scatter vector for the element
    sctrx=sctr;           % x scatter vector
    sctry=sctrx+numnode;  % y scatter vector
    for q=1:size(W,1)
        pt=Q(q,:);
        wt=W(q);
        [N,dNdxi]=lagrange_basis(edgeElemType,pt);  % element shape functions
        J0=dNdxi'*node(sctr,:);
        detJ0=norm(J0);
        % element Jacobian
        % determiniat of jacobian
        yPt=N'*node(sctr,2);
        fyPt=-P*(c^2-yPt^2)/(2*I0);
        %f(sctry)=f(sctry)+N*fyPt*detJ0*wt;  %tor
        f(sctrx)=f(sctrx)+N*1000*detJ0*wt;
    end % of quadrature loop
end  % of element loop


%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'   COMPUTING STIFFNESS MATRIX'])
switch elemType  % define quadrature rule
    case 'Q9'
        [W,Q]=quadrature( 4, 'GAUSS', 2 ); % 4x4 Gaussian quadrature
    case 'Q4'
        [W,Q]=quadrature( 2, 'GAUSS', 2 ); % 2x2 Gaussian quadrature
    otherwise
        [W,Q]=quadrature( 1, 'TRIANGULAR', 2 ); % 1 point triangural quadrature
end
for e=1:numelem                          % start of element loop
    sctr=element(e,:);           % element scatter vector
    sctrB=[ sctr sctr+numnode ]; % vector that scatters a B matrix
    nn=length(sctr);
    for q=1:size(W,1)
        pt=Q(q,:);
        wt=W(q);
        [N,dNdxi]=lagrange_basis(elemType,pt); % element shape functions
        J0=node(sctr,:)'*dNdxi;                % element Jacobian matrix
        invJ0=inv(J0);
        dNdx=dNdxi*invJ0;
        
        B(1,1:nn)       = dNdx(:,1)';
        B(2,nn+1:2*nn)  = dNdx(:,2)';
        B(3,1:nn)       = dNdx(:,2)';
        B(3,nn+1:2*nn)  = dNdx(:,1)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
        K(sctrB,sctrB)=K(sctrB,sctrB)+B'*C*B*W(q)*det(J0);
    end  % of quadrature loop
end    % of element loop
%%%%%%%%%%%%%%%%%%% END OF STIFFNESS MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%%%%

% APPLY ESSENTIAL BOUNDARY CONDITIONS
disp([num2str(toc),'   APPLYING BOUNDARY CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average size of an element in K
% used to keep the conditioning of the K matrix
udofs=fixedNode;           % global indecies of the fixed x displacements
vdofs=fixedNode+numnode;   % global indecies of the fixed y displacements
f=f-K(:,udofs)*uFixed;  % modify the force vector
f=f-K(:,vdofs)*vFixed;
K(udofs,:)=0;
K(vdofs,:)=0;
K(:,udofs)=0;
K(:,vdofs)=0;
K(udofs,udofs)=bcwt*speye(length(udofs)); % put ones*bcwt on the diagonal
K(vdofs,vdofs)=bcwt*speye(length(vdofs));
f(udofs)=bcwt*speye(length(udofs))*uFixed;
f(vdofs)=bcwt*speye(length(udofs))*vFixed;
% SOLVE SYSTEM
disp([num2str(toc),'   SOLVING SYSTEM'])
U=K\f;
%******************************************************************************
%*** POST - PROCESSING *** %***************************************************

ptId  = intersect (find(node(:,1)==48),find(node(:,2)==0));
numUy = U(ptId+numnode)

% Here we plot the stresses and displacements of the solution. As with the
% mesh generation section we don?t go into too much detail - use help
% ?function name? to get more details.
disp([num2str(toc),'   POST-PROCESSING'])
dispNorm=L/max(sqrt(U(xs).^2+U(ys).^2));
scaleFact=0.1*dispNorm;
fn=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DEFORMED DISPLACEMENT PLOT
Ux = U(xs);
Uy = U(ys);

figure(fn)
clf
plot_field(node+scaleFact*[Ux Uy],element,elemType,U(xs));
hold on
plot_mesh(node+scaleFact*[Ux Uy],element,elemType,'g.-');
%plot_mesh(node,element,elemType,'w--');
colorbar
fn=fn+1;
title('DEFORMED DISPLACEMENT IN Y-DIRECTION')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE STRESS
stress=zeros(numelem,size(element,2),3);
switch elemType  % define quadrature rule
    case 'Q9'
        stressPoints=[-1 -1;1 -1;1 1;-1 1;0 -1;1 0;0 1;-1 0;0 0 ];
    case 'Q4'
        stressPoints=[-1 -1;1 -1;1 1;-1 1];
    otherwise
        stressPoints=[0 0;1 0;0 1];
end

for e=1:numelem
    sctr=element(e,:);
    sctrB=[sctr sctr+numnode];
    nn=length(sctr);
    for q=1:nn
        pt=stressPoints(q,:);
        [N,dNdxi]=lagrange_basis(elemType,pt);
        J0=node(sctr,:)'*dNdxi;
        invJ0=inv(J0);
        dNdx=dNdxi*invJ0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE B MATRIX
        B=zeros(3,2*nn);
        B(1,1:nn)       = dNdx(:,1)';
        B(2,nn+1:2*nn)  = dNdx(:,2)';
        B(3,1:nn)       = dNdx(:,2)';
        B(3,nn+1:2*nn)  = dNdx(:,1)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE ELEMENT STRAIN AND STRESS AT STRESS POINT
        strain=B*U(sctrB);
        stress(e,q,:)=C*strain;
    end
end   % of element loop
stressComp=1;
figure(fn)
clf
plot_field(node+scaleFact*[U(xs) U(ys)],element,elemType,stress(:,:,stressComp));
hold on
%plot_mesh(node+scaleFact*[U(xs) U(ys)],element,elemType,'g.-');
%plot_mesh(node,element,elemType,'w--');
colorbar
fn=fn+1;
title('DEFORMED STRESS PLOT, BENDING COMPONENT')
%print(fn,?-djpeg90?,[?beam_?,elemType,?_sigma?,num2str(stressComp),?.jpg?])

computeFemEnergyNorm

% Timoshenko beam problem: Cantilever beam subjected to a parabolic force
clear all
close all

addpath ../nurbs-util/


%%%%%%%%%%%%%%%%%%
% INPUT DATA
%%%%%%%%%%%%%%%%%

%Material properties
E = 3e4; 
nu = 0.3;

%dimensions of the beam
L = 48;  %length of the beam
W = 12;  %width of the beam

numix = 12; % number of elements along the length
numiy = 3; % number of elements along the width
%NOTE: For best results, numix, numiy should be proportional to L,W

%order of Nurbs in each direction
p = 2;  %p=q=1 implies Q4 FEM
q = 2;

%parameter P in Timoshenko beam problem
P = 10;

%%%%%%%%%%%%%%%%%%%%%%%%
% END INPUT DATA
%%%%%%%%%%%%%%%%%%%%%%%%


tic
disp('Preprocessing...')

%construct the knot vectors, repeating the endpoints p+1 times
knotu = [zeros(1,p), linspace(0,1,numix+1), ones(1,p)];
knotv = [zeros(1,q), linspace(0,1,numiy+1), ones(1,q)];

%the number of control points in the u and v directions
lenu = length(knotu)-p-1;  %number of basis functions in the u direction
lenv = length(knotv)-q-1;  %number of basis functions in the v direction
numnodes = lenu*lenv;      %number of control points

ngaussx = p+1; %number of Gauss points used in the x-direction
ngaussy = q+1; %number of Gauss points used in the y-direction
ngaussedge = max(p+1,q+1); %number of Gauss points used for edge integration

dim = 2; % we are on a two-dimensional plate

%intialize the control net
coordinates = zeros(numnodes, dim+1);  %allocate one more column for the weights
coord_ij = zeros(numnodes, dim); %coordinate # -> i,j array
b_net = zeros(lenu, lenv, dim+1); %same array in the B_{ij} format

for j=1:lenv %for each node in the y direction
    
    %compute the control points in the x-direction using Greville Abscisae
    coordy = W.*sum(knotv(j+1:j+q))./q;
      
    for i=1:lenu % for each node in the x direction
        index = (j-1)*lenu + i; %the index of the coordinate array
        coordx = L.*sum(knotu(i+1:i+p))./p;
        coordinates(index,:) = [coordx, coordy, 1]; %put the (i,j) node in the coordinate array with weight 1      
        b_net(i,j,:) = coordinates(index,:);        
        coord_ij(index,:) = [i, j];
    end %for i
end %for j

%initialize the element integration matrix
luku = length(unique(knotu))-1;
lukv = length(unique(knotv))-1;
element_int = zeros(luku*lukv,4); %corners in parameter space of each element (knot-span) [(u_i, u_{i+1}]x[v_j, v_{j+1}]
nument = (p+1)*(q+1); %number of nodes for each element
element_nod = zeros(luku*lukv,nument); %element -> node connectivity

%loop through each element and compute element_int, element_nod
elementcounter = 0;
for j=1:length(knotv)-1
    for i=1:length(knotu)-1
        if (knotu(i+1)>knotu(i) && knotv(j+1) >knotv(j))  %the knotspan has non-zero area
%            [i, j]
            elementcounter = elementcounter + 1;
            element_int(elementcounter,:)= [knotu(i), knotv(j), knotu(i+1), knotv(j+1)];
            tcount = 0;
            currow = zeros(1, nument);
            %now we add the nodes from i-p...i in the u direction and
            %j-q...j in the v direction
            for t2=j-q:j
                for t1 = i-p:i                                        
                    tcount = tcount + 1;                    
                    currow(tcount) = t1+(t2-1)*lenu;
                end
            end
            element_nod(elementcounter,:)=currow;
        end
    end
end


%find out the nodes of functions that have support on each segment of edge
%format: vertex1 vertex2 element# node1 node2 ...

bottomedge = zeros(numix, 5);
%bottomcount = 0;

%for each edge, pick up the coordinates in parameter space from element int
%orientation: 1-bottom, 2-right, 3-top, 4-left
for t=1:numix    
    elt = t;    
    %format: nodes = [ximin, ximax, eta0];
    nodes = [element_int(elt,1), element_int(elt,3), element_int(elt,2)];
    ort = 1;
    bottomedge(t, :) = [elt, nodes, ort];
        
end

topedge = zeros(numix, 5);
for t=1:numix        
    elt = elementcounter-numix+t;
    %format: nodes = [ximin, ximax, eta0];
    nodes = [element_int(elt,1), element_int(elt,3), element_int(elt,4)];
    ort = 3;    
    topedge(t,:) = [elt, nodes, ort];
end

leftedge = zeros(numiy, 5);
for t=1:numiy
    vertices = [(t-1)*lenu + 1, t*lenu+1];    
    %format: nodes = [etamin, etamax, xi0];
    elt = (t-1)*numix+1;
    nodes = [element_int(elt,2), element_int(elt,4), element_int(elt,1)];    
    ort = 4;
    leftedge(t,:) = [elt, nodes, ort];
end


rightedge = zeros(numiy, 5);
for t=1:numiy    
    %format: nodes = [etamin, etamax, xi0];
    elt = t*numix;
    nodes = [element_int(elt,2), element_int(elt,4), element_int(elt,3)];   
    ort = 2;
    rightedge(t,:) = [elt, nodes, ort];
end

%pick wich edges correspond to neumann and dirichlet conditions
neumann = rightedge;
dirichlet = leftedge;
dirichletx = leftedge;
dirichlety = leftedge;

globnum = dim*lenu*lenv; %number of global DOF

stiff = sparse(globnum, globnum);
rhs = zeros(globnum,1);

%assembly
toc
disp('Assembling stiffness matrix and right hand side...')

%calculate the elasticity matrix
C=E/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];  %plane stress

for i=1:elementcounter
    
    ximin = element_int(i,1);
    etamin = element_int(i,2);
    ximax = element_int(i,3);
    etamax = element_int(i,4);            

    scalefac = (ximax - ximin)*(etamax - etamin)/4;
    [gpx,gwx]=genGP_GW(ngaussx);
    [gpy,gwy]=genGP_GW(ngaussy);           
  
    localstiff = zeros(dim*nument, dim*nument);
    localrhs = zeros(dim*nument, 1);
    scrtx = zeros(dim*nument, 1);
    %for each quadrature point evaluate the shape functions and calculuate
    %local stiffness matrix and local rhs vector
    for ii=1:ngaussx
        for jj=1:ngaussy
             %compute the value of the shape function and gradient                          
            [R, dRdx, ~, J] = nurbshaped(i, gpx(ii), gpy(jj), knotu, knotv, b_net, p, q, lenu, lenv, element_nod, coord_ij);                        
            B = zeros(dim*nument, 3);
            %loop over each shape function and build the B matrix
            for j=1:nument
                globnum = element_nod(i,j);                                
                cR = R(j);
                cdRdx = dRdx(1,j);
                cdRdy = dRdx(2,j);
                
                scrtx(2*j-1) = 2*globnum-1;
                scrtx(2*j) = 2*globnum;
                B(2*j-1, :) = [cdRdx, 0, cdRdy];
                B(2*j, :) = [0, cdRdy, cdRdx];
               
                %there is zero volume force for beam problem, so we don't
                %do anything (localrhs remains zero)
            end           
            localstiff = localstiff + B*C*B'.*J.*scalefac.*gwx(ii).*gwy(jj);            
        end
    end
    stiff(scrtx,scrtx) = stiff(scrtx, scrtx) + localstiff;
    rhs(scrtx) = rhs(scrtx) + localrhs;  
end


%impose Neumann conditions
toc
disp('Imposing Neumann conditions...')
for i=1:size(neumann,1)
    
    elementind = neumann(i,1);
    
    curnod = element_nod(elementind, nument);
    ni = coord_ij(curnod, 1);
    nj = coord_ij(curnod, 2);
    
    corient = neumann(i, 5);
    
    if (corient == 1) || (corient == 3)
        scalefac = (knotu(ni+1) - knotu(ni))/2;
    else
        scalefac = (knotv(nj+1) - knotv(nj))/2;
    end
    [gp, gw] = genGP_GW(ngaussedge);
    localrhsed = zeros(dim*nument, 1);
    scrtx = zeros(dim*nument, 1);

    %loop over Gauss points and compute the integral
    for igauss = 1:ngaussedge
        [R, coord, normal, J] = nurbedge(gp(igauss), ni, nj, knotu, knotv, b_net, p, q, lenu, lenv, corient);
      
        gwt = gw(igauss)*scalefac;       
        [~, stress] = gbeam(coord(1), coord(2), P, E, nu, W, L);        
        
        for j=1:nument
            globnum = element_nod(elementind,j);                                                 
            taux = normal(1)*stress(1) + normal(2)*stress(3);
            tauy = normal(1)*stress(3) + normal(2)*stress(2);
                       
            locnum = nument+1-j;
            cR = R(locnum);
            
            scrtx(2*j-1) = 2*globnum-1;
            scrtx(2*j) = 2*globnum;
                     
            localrhsed(2*j-1) = localrhsed(2*j-1) + taux.*cR.*J.*gwt;
            localrhsed(2*j) = localrhsed(2*j) + tauy.*cR.*J.*gwt;
        end        
    end
    rhs(scrtx)=rhs(scrtx)+localrhsed;
end


toc
disp('Imposing Dirichlet boundary conditions')

%old method of imposing Dirichlet boundary condtions on the left edge
%uses interpolation at Control Points

DirichletNodes = 1:lenu:lenu*(lenv-1)+1;
lendn = length(DirichletNodes);
globaldir = zeros(dim*lendn,1);

%create an array of Dirichlet nodes in parameter space

DirichletParmNodes = linspace(0,1,lendn);
DirichletOrient = 4; %left edge

%set-up local-global connectivity matrix for interpolation
if DirichletOrient == 4 
    loc2globint = zeros(nument, length(knotv)-1);
    for i=1:length(knotv)-1
        loc2globint(nument:-(p+1):p+1,i)=(i-p):i;
    end
end
        
InterpolationMatrix = zeros(lendn,lendn);
InterpolationRHSX = zeros(lendn, 1);
InterpolationRHSY = zeros(lendn, 1);

for i=1:lendn  %loop over the sample points
    curnode = DirichletParmNodes(i);
    if DirichletOrient == 4
        knot = knotv;
        ni = p+1;
        for j=1:length(knot)-1
            
            if (knot(j)<knot(j+1)) && (knot(j)<=curnode) && (curnode <= knot(j+1))
                nj = j;
                %convert to master coordinates on interval [-1, 1]
                coord = (2*curnode-knot(j)-knot(j+1))/(knot(j+1)-knot(j));
                [shb, phyc, normal, detjb] = nurbedge(coord, ni, nj, knotu, knotv, b_net, p, q, lenu, lenv, DirichletOrient);
                
                for k = p+1:p+1:nument                    
                    globnum = loc2globint(k, nj);
                    InterpolationMatrix(i,globnum) = shb(k);
                end

                NodeEval = gbeam(phyc(1),phyc(2),P, E, nu, W, L);       
                InterpolationRHSX(i) = NodeEval(1);
                InterpolationRHSY(i) = NodeEval(2);
            end
        end
    end
end
                
InterpolationValuesX = InterpolationMatrix\InterpolationRHSX;
InterpolationValuesY = InterpolationMatrix\InterpolationRHSY;

%impose Dirichlet conditions for x-displacement
[A,b]=feaplyc2(stiff,rhs,2*DirichletNodes-1,InterpolationValuesX);
%impose Dirichlet conditions for y-displacement
[A,b]=feaplyc2(A,b,2*DirichletNodes,InterpolationValuesY);

toc
disp('Solving the linear system...')
conditionest = condest(A)
%condition = cond(full(A))
% Calculating the solution
x = A\b;
u = x;
maxresidual = max(abs(A*x-b))
%energynorm = sqrt(x'*A*x)
%EU = 0.5*u'*stiff*u;
toc
disp('Visualizing the result...')
size(A)
visualize

figure
plotBoundary

toc
disp('Calculating the error...')
l2err
toc



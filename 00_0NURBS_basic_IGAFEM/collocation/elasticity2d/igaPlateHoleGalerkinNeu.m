%solves the "Plate with a hole" problem with Galerkin method
%Neumann boundary condtion + symmetry condtions

%the plate ia centered at the origin, the quarter-plate modeled lies in
%the IInd quadrant.


close all
clear all

addpath ../../nurbs-geopdes/inst  %to use functions from Octave NURBS package

L = 4;  %length of the plate
rad = 1;  %radius of the hole
tx = 10; %force applied in the x direction
ty = 0;  %force applied in the y direction


[Emod, nu] = materialenu(0.5,0.5);
C=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];


for p=5:5 %degree of Nurbs in each direction
    q=p;
    for numrefinements = 4:4
        tic
        disp('Preprocessing...')
        
        dim = 2; % we are on a two-dimensional plate
        
        %define geometry on coarsest mesh
        knotu = [0, 0, 0, 0.5, 1, 1, 1];
        knotv = [0, 0, 0, 1, 1, 1];
        
        weights = [1,(1+1/sqrt(2))/2,(1+1/sqrt(2))/2,1, 1,1,1,1, 1,1,1,1]';
        controlPoints = [-1, 0; -1, sqrt(2)-1; 1-sqrt(2), 1; 0, 1; -2.5, 0; -2.5, 0.75; -0.75, 2.5; 0, 2.5; -4, 0; -4, 4; -4, 4; 0, 4];
        
        %the number of control points in the u and v directions
        lenu = length(knotu)-3;  %number of basis functions in the  u direction
        lenv = length(knotv)-3;  %number of basis functions in the v direction
        coordinates = [controlPoints, weights];
        
        %convert control points/weights to coefs format in NURBS toolbox
        coefs = zeros(4,lenu,lenv);
        for i=1:lenu
            for j=1:lenv
                index = (j-1)*lenu+i;
                xcoord = coordinates(index,1);
                ycoord = coordinates(index,2);
                wght = coordinates(index,3);
                coefs(1,i,j) = xcoord*wght;
                coefs(2,i,j) = ycoord*wght;
                coefs(4,i,j) = wght;
            end
        end
        knots =  {knotu, knotv};
        nurbs = nrbmak(coefs, knots); %generates the coarsest nurbs object, p=2, q=2
        
        nurbsr = nrbdegelev(nurbs, [p,q]-(nurbs.order-1)); %elevate to given p,q
        
        %convert NURBS toolbox format back to control points/weights
        knotu = nurbsr.knots{1};
        knotv = nurbsr.knots{2};
        
        %calculate knots to be inserted during h-refinement
        uknotu = unique(knotu);
        uknotv = unique(knotv);
        for i=1:numrefinements
            %take the average of the knots and store them in insvecu,
            %insvecv
            
            insvec1 = uknotu(1:end-1);
            insvec2 = uknotu(2:end);
            insvecu = (insvec1+insvec2)/2;
            
            insvec1 = uknotv(1:end-1);
            insvec2 = uknotv(2:end);
            insvecv = (insvec1+insvec2)/2;
            
            
            uknotu = sort([uknotu, insvecu]);
            uknotv = sort([uknotv, insvecv]);
            
        end
        
        insvecu = setdiff(uknotu,[0,1/2,1]);
        insvecv = uknotv(2:end-1);
        
        nurbsr = nrbkntins(nurbsr, {insvecu, insvecv}); %insert the new knot points
        knotu = nurbsr.knots{1};
        knotv = nurbsr.knots{2};
        
        %the number of control points in the u and v directions
        lenu = length(knotu)-p-1;  %number of basis functions in the u direction
        lenv = length(knotv)-q-1;  %number of basis functions in the v direction
        
        controlPoints = zeros(lenu*lenv,2);
        weights = zeros(lenu*lenv,1);
        
        for i=1:lenu
            for j=1:lenv
                index = (j-1)*lenu+i;
                xcoord = nurbsr.coefs(1,i,j)/nurbsr.coefs(4,i,j);
                ycoord = nurbsr.coefs(2,i,j)/nurbsr.coefs(4,i,j);
                wght = nurbsr.coefs(4,i,j);
                controlPoints(index,:) = [xcoord, ycoord];
                weights(index) = wght;
            end
        end
        
        coordinates = [controlPoints, weights];
        
        %the number of control points in the u and v directions
        lenu = length(knotu)-p-1;  %number of basis functions in the u direction
        lenv = length(knotv)-q-1;  %number of basis functions in the v direction
        numnodes = lenu*lenv;      %number of control points
        numix = length(unique(knotu))-1;
        numiy = length(unique(knotv))-1;
        
        if p==2
            ngaussx = p+2; %number of Gauss points used in the x-direction
            ngaussy = q+2; %number of Gauss points used in the y-direction
        else
            ngaussx = p+1;
            ngaussy = p+1;
        end
        
        
        %ngaussedge = max(ngaussx,ngaussy); %number of Gauss points used for edge integration (Dirichlet L2 projection or Neumann conditions)
        ngaussedge = 8;
        
        %intialize the control net
        %coordinates = zeros(numnodes, dim+1);  %allocate one more column for the weights
        coord_ij = zeros(numnodes, dim); %coordinate # -> i,j array
        b_net = zeros(lenu, lenv, dim+1); %same array in the B_{ij} format
        
        for j=1:lenv %for each node in the y direction
            for i=1:lenu % for each node in the x direction
                index = (j-1)*lenu + i; %the index of the coordinate array
                b_net(i,j,:) = coordinates(index,:);
                coord_ij(index,:) = [i, j];
            end %for i
        end %for j
        
        %initialize the element integration matrix
        %format of each row: <#of integration subdomains> <type of subdomain>
        %[x1][y1][x2][y2][<type of subdomain>][x1][y1]...
        luku = length(unique(knotu))-1;
        lukv = length(unique(knotv))-1;
        
        nument = (p+1)*(q+1); %number of nodes for each element
        element_int = zeros(luku*lukv,4);
        element_nod = zeros(luku*lukv,nument);
        tolarea = 1e-8;
        
        
        elementcounter = 0;
        for j=1:length(knotv)-1
            for i=1:length(knotu)-1
                if (abs(knotu(i+1)-knotu(i))>tolarea) && (abs(knotv(j+1)-knotv(j))>tolarea)  %the knotspan has non-zero area
                    %            [i, j]
                    elementcounter = elementcounter + 1;
                    element_int(elementcounter,:) = [knotu(i), knotv(j), knotu(i+1), knotv(j+1)];
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
        
        
        %for each edge, pick up the coordinates in parameter space from element int
        %orientation: 1-bottom, 2-right, 3-top, 4-left
        bottomedge = zeros(numix, 5);
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
        
        %pick wich edges correspond to dirichlet conditions
        
        dirichletx = rightedge;
        dirichlety = leftedge;
        neumann = topedge;
        
        toc
        disp('Pre-computing B-Splines...')
        
        %loop throught the knot spans in each direction and pre-compute
        %Bsplines and connectivity arrays
        
        [gpx,gwx]=genGP_GW(ngaussx);
        [gpy,gwy]=genGP_GW(ngaussy);
        
        deriv_order = 1; %we only need 1st derivatives for Galerkin method
        [ pt_index_u, M_arr_u ] = makeIndexBspline( knotu, p, deriv_order, gpx );
        [ pt_index_v, M_arr_v ] = makeIndexBspline( knotv, q, deriv_order, gpy );
        
        
        toc
        disp('Assembling stiffness matrix and right hand side...')
        
        %index for storing entries of the LHS
        %element_block_size = elementcounter;
        element_block_size = numix;
        II = zeros(1, dim^2*element_block_size*nument^2);
        JJ = zeros(1, dim^2*element_block_size*nument^2);
        indexcounter = 0;
        S = zeros(1, dim^2*element_block_size*nument^2);
        stiff = sparse(dim*numnodes, dim*numnodes);
        rhs = zeros(dim*numnodes,1);
        %assembly
        
        
        %for each integration element calculate the value of the basis
        %functions
        for i=1:elementcounter
            ximin = element_int(i,1);
            etamin = element_int(i,2);
            ximax = element_int(i,3);
            etamax = element_int(i,4);
            
            scalefac = (ximax - ximin)*(etamax - etamin)/4;
            
            %find the knot span index in each direction
            ni = coord_ij(element_nod(i,end),1);
            nj = coord_ij(element_nod(i,end),2);
            
            localstiff = zeros(2*nument, 2*nument); %local stiffness
            
            localrhs = zeros(nument, 1);
            scrtx = element_nod(i,:);
            
            %calculate the weights and control points corresponding to the
            %current element
            
            wgts = reshape(b_net(ni-p:ni, nj-q:nj, dim+1), nument, 1);
            cpts = reshape(b_net(ni-p:ni, nj-q:nj, 1:2), nument, dim);
            
            %for each shape function in shapelist compute the matrices B, E
            for ii=1:ngaussx
                for jj=1:ngaussy
                    
                    glob_index_u = pt_index_u(ni, ii);
                    glob_index_v = pt_index_v(nj, jj);
                    
                    M = M_arr_u(:,:,glob_index_u);
                    P = M_arr_v(:,:,glob_index_v);
                    B = zeros(2*nument,3);
                    
                    [R, dRdx] = nurbshape2d4(M,P,p,q,wgts);
                    
                    %calculate the coordinates in the physical space
                    coord = R*cpts;
                    
                    %calculate the Jacobian of the transformation
                    dxdxi = dRdx*cpts;
                    J = det(dxdxi);
                    
                    %calculate the derivatives with respect to physical
                    %space
                    dRdx = dxdxi\dRdx;
                    dRdx = dRdx';
                    
                    B(1:2:2*nument-1,1) = dRdx(:,1);
                    B(2:2:2*nument,2) = dRdx(:,2);
                    B(1:2:2*nument-1,3) = dRdx(:,2);
                    B(2:2:2*nument,3) = dRdx(:,1);
                    
                    %TODO: implement non-zero volume force
                    localstiff = localstiff + B * C * B' * J * scalefac * gwx(ii).*gwy(jj);
                    
                end
            end
            
            %double the entries in scrtx
            dscrtx = reshape([2*scrtx-1; 2*scrtx],1,2*nument);
            if (indexcounter+dim^2*nument^2)>element_block_size*nument^2*dim^2
                stiff = stiff + sparse(II,JJ,S,dim*numnodes,dim*numnodes);
                indexcounter = 0;
                % disp(sprintf('Element %d of %d',i, elementcounter))
                % toc
            end
            
            II(indexcounter+1:indexcounter+dim^2*nument^2) = repmat(dscrtx,1,dim*nument);
            JJ(indexcounter+1:indexcounter+dim^2*nument^2) = reshape(repmat(dscrtx',1,dim*nument)',1,dim^2*nument^2);
            S(indexcounter+1:indexcounter+dim^2*nument^2) = reshape(localstiff,1,dim^2*nument^2);
            indexcounter = indexcounter + dim^2*nument^2;
            
            rhs(scrtx) = rhs(scrtx) + localrhs;
            
        end
        
        stiff = stiff + sparse(II,JJ,S,dim*numnodes,dim*numnodes);
        clear II JJ S
        
        toc
        %impose Neumann conditions
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
            localrhsed = zeros(2*nument, 1);
            scrt = zeros(2*nument, 1);
            
            %loop over Gauss points and compute the integral
            for igauss = 1:ngaussedge
                [R, coord, normal, J] = nurbedge(gp(igauss), ni, nj, knotu, knotv, b_net, p, q, corient);
                gwt = gw(igauss)*scalefac;
                stress = ghole(coord(1), coord(2), rad, tx);
                
                for t=1:nument
                    
                    taux = normal(1)*stress(1) + normal(2)*stress(3);
                    tauy = normal(1)*stress(3) + normal(2)*stress(2);
                    
                    cR = R(end+1-t);
                    scrt(2*t-1:2*t) = [2*element_nod(elementind,t)-1, 2*element_nod(elementind,t)];
                    
                    localrhsed(2*t-1) = localrhsed(2*t-1) + taux.*cR.*J.*gwt;
                    localrhsed(2*t) = localrhsed(2*t) + tauy.*cR.*J.*gwt;
                end
                %pause
            end
            rhs(scrt)=rhs(scrt)+localrhsed;
            
        end
        
        
        toc
        disp('Imposing Dirichlet boundary conditions...')
        
        [dirichlet_nodesx, dirvalsx] = DirichletProj(dirichletx, 1, element_nod, numnodes, knotu, knotv, ngaussedge, nument, b_net, p, q, lenu, lenv, coord_ij, coordinates, Emod, nu, tx, rad);
        [dirichlet_nodesy, dirvalsy] = DirichletProj(dirichlety, 2, element_nod, numnodes, knotu, knotv, ngaussedge, nument, b_net, p, q, lenu, lenv, coord_ij, coordinates, Emod, nu, tx, rad);
        
        %impose the boundary conditions directly on the linear system
        [stiff,rhs]=feaplyc2(stiff,rhs,2*dirichlet_nodesx-1,dirvalsx);
        [stiff,rhs]=feaplyc2(stiff,rhs,2*dirichlet_nodesy,dirvalsy);
        
        
        toc
        disp('Solving the linear system...')
        
        % Calculating the solution
        x = stiff\rhs;
        
        toc
        disp('Calculating the error...')
        [l2errrel, relenergynorm] = h1err(x, element_nod, ngaussx, ngaussy, elementcounter, element_int, coord_ij, knotu, knotv, b_net, p, q, C, rad, Emod, nu, tx);
        fprintf('p= %d num_ref = %d dof = %d l2err = %1.15f energy_err= %1.15f \n', p, numrefinements, length(x), l2errrel, relenergynorm)
        
        toc
        
    end
end

disp('Visualizing the result...')
visualize
toc

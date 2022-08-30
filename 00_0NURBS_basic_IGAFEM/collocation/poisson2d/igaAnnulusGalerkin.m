% script igaAnnulusGalerkin.m
% Solve -Delta(U) + alpha*u = f on a quarter annulus domain with homogeneous Dirichlet boundary
% Galerkin method
% INPUT data: p - polynomial order
%             numiset - set that contains the number of elements in each
%             direction


addpath ../../nurbs-geopdes/inst
close all
clear all

numiset = [10:10:20]; % number of intervals set
[~, ~, ~, ~, alpha] = exact_sol(0,0); %solve -Delta u + alpha*y = f.

tolcomp = 1e-10; %tolerance for detecting non-empty knot spans

for p=4:4
    q = p;
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [1;0;0];
    coefs(1:3,1,2) = [sqrt(2)/2; sqrt(2)/2; 0];
    coefs(1:3,1,3) = [0;1;0];
    coefs(1:3,2,1) = [4;0;0];
    coefs(1:3,2,2) = [2*sqrt(2);2*sqrt(2);0];
    coefs(1:3,2,3) = [0;4;0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = sqrt(2)/2;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = sqrt(2)/2;
    coefs(4,2,3) = 1;
    
    knots = {[0 0 1 1], [0 0 0 1 1 1]};
    nurbs = nrbmak(coefs, knots); %generates the coarsest nurbs object, p=1, q=2
    nurbs = nrbdegelev(nurbs, [p,q]-(nurbs.order-1)); %elevate to given p,q
    
    for numix = numiset % number of elements along the length
        
        numiy = numix; %set same number of elements along the width
        tic
        disp('Pre-processing...')
        
        %refine by knot insertion
        knot1 = linspace(0,1,numix+1);
        knot2 = linspace(0,1,numiy+1);
        
        nurbsr = nrbkntins(nurbs, {knot1(2:end-1) knot2(2:end-1)});
        knotu = nurbsr.knots{1};
        knotv = nurbsr.knots{2};
        
        
        %the number of control points in the u and v directions
        lenu = length(knotu)-p-1;  %number of basis functions in the u direction
        lenv = length(knotv)-q-1;  %number of basis functions in the v direction
        numnodes = lenu*lenv;      %number of control points
        
        dim = 2; % we are on a two-dimensional world
        
        if p==3
            ngaussx = p+2;
            ngaussy = q+2;
        else
            ngaussx = p+1;
            ngaussy = q+1;
        end
        
        ngaussedge = max(ngaussx,ngaussy);
        
        
        %intialize the control net
        coordinates = zeros(numnodes, dim+1);  %allocate one more column for the weights
        coord_ij = zeros(numnodes, dim); %coordinate # -> i,j array
        b_net = zeros(lenu, lenv, dim+1); %same array in the B_{ij} format
        
        for j=1:lenv %for each node in the y direction
            for i=1:lenu % for each node in the x direction
                index = (j-1)*lenu + i; %the index of the coordinate array
                coordinates(index,:) = [nurbsr.coefs(1,i,j)./nurbsr.coefs(4,i,j), nurbsr.coefs(2,i,j)./nurbsr.coefs(4,i,j), nurbsr.coefs(4,i,j)]; %put the (i,j) node in the coordinate array
                b_net(i,j,:) = coordinates(index,:);
                coord_ij(index,:) = [i, j];
            end %for i
        end %for j
        
        
        luku = length(unique(knotu))-1;
        lukv = length(unique(knotv))-1;
        
        nument = (p+1)*(q+1); %number of nodes for each element
        element_nod = zeros(luku*lukv,nument);
        element_int = zeros(luku*lukv,4);
        
        %loop through each element and compute element_int, element_nod
        
        elementcounter = 0;
        for j=1:length(knotv)-1
            for i=1:length(knotu)-1
                if (abs(knotu(i+1)-knotu(i))>tolcomp) && (abs(knotv(j+1)-knotv(j))>tolcomp)  %the knotspan has non-zero area
                    % [i, j]
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
        
        neumann = [];
        dirichlet = [leftedge;rightedge;bottomedge;topedge];
        
        toc
        disp('Pre-computing B-Splines...')
        
        %loop through the knot spans in each direction and pre-compute
        %Bsplines and connectivity arrays
        
        [gpx,gwx]=genGP_GW(ngaussx);
        [gpy,gwy]=genGP_GW(ngaussy);
        
        deriv_order = 1; %we only need 1st derivatives for Galerkin method
        [ pt_index_u, M_arr_u ] = makeIndexBspline( knotu, p, deriv_order, gpx );
        [ pt_index_v, M_arr_v ] = makeIndexBspline( knotv, q, deriv_order, gpy );
        %assembly
        
        toc
        disp('Assembling linear system...')
        
        %index for storing entries of the LHS
        element_block_size = elementcounter;
        II = zeros(1, element_block_size*nument^2);
        JJ = zeros(1, element_block_size*nument^2);
        indexcounter = 0;
        S = zeros(1, element_block_size*nument^2);
        %stiff = sparse(numnodes, numnodes);
        rhs = zeros(numnodes, 1);
        
        for i=1:elementcounter
            
            %for each integration element calculate the value of the basis
            %functions
            
            ximin = element_int(i,1);
            etamin = element_int(i,2);
            ximax = element_int(i,3);
            etamax = element_int(i,4);
            
            scalefac = (ximax - ximin)*(etamax - etamin)/4;
            
            %find the knot span index in each direction
            ni = coord_ij(element_nod(i,end),1);
            nj = coord_ij(element_nod(i,end),2);
            
            localstiff = zeros(nument, nument); %local stiffness and mass matrix
            
            localrhs = zeros(nument, 1);
            scrtx = element_nod(i,:);
            
            %calculate the weights and control points corresponding to the
            %current element
            
            wgts = reshape(b_net(ni-p:ni, nj-q:nj, dim+1), nument, 1);
            cpts = reshape(b_net(ni-p:ni, nj-q:nj, 1:2), nument, dim);
            
            %for each shape function in shapelist compute the matrices B, E
            for ii=1:ngaussx
                for jj=1:ngaussy
                    %[R_old, dRdx_old, ~, coord_old, J_old, dxdxi_old] = nurbshaped(i, gpx(ii), gpy(jj), knotu, knotv, b_net, p, q, lenu, lenv, element_nod, coord_ij);
                    
                    glob_index_u = pt_index_u(ni, ii);
                    glob_index_v = pt_index_v(nj, jj);
                    
                    M = M_arr_u(:,:,glob_index_u);
                    P = M_arr_v(:,:,glob_index_v);
                    
                    [R, dRdx] = nurbshape2d4(M,P,p,q,wgts);
                    
                    %calculate the coordinates in the physical space
                    coord = R*cpts;
                    
                    %calculate the Jacobian of the transformation
                    dxdxi = dRdx*cpts;
                    J = det(dxdxi);
                    
                    %calculate the derivatives with respect to physical
                    %space
                    dRdx = dxdxi\dRdx;
                    
                    [~, ~, ~, f] = exact_sol(coord(1),coord(2));
                    localrhs = localrhs + f.*R'.*J.*scalefac.*gwx(ii).*gwy(jj);
                    localstiff = localstiff + (dRdx'*dRdx+alpha*(R'*R)).*J.*scalefac.*gwx(ii).*gwy(jj);
                    
                end
            end
            
            II(indexcounter+1:indexcounter+nument^2) = repmat(scrtx,1,nument);
            JJ(indexcounter+1:indexcounter+nument^2) = reshape(repmat(scrtx',1,nument)',1,nument^2);
            S(indexcounter+1:indexcounter+nument^2) = reshape(localstiff,1,nument^2);
            indexcounter = indexcounter + nument^2;
            
            
            rhs(scrtx) = rhs(scrtx) + localrhs;
            
            
        end
        
        stiff = sparse(II,JJ,S,numnodes, numnodes);
        clear II JJ S
        toc
        
        disp('Imposing Dirichlet boundary conditions...')
        [bcdof] = DirichletProj(dirichlet, element_nod, numnodes, knotu, knotv, ngaussedge, nument, b_net, p, q, lenu, lenv, coord_ij, coordinates);
        bcval = zeros(size(bcdof));
        [stiff,rhs] = feaplyc2_poisson(stiff,rhs,bcdof);
        toc
        
        disp('Solving the linear system...')
        
        sol0 = stiff\rhs;
        
        
        [bcdofs, ind] = sort(bcdof);
        bcvals = bcval(ind);
        for i=1:length(bcdof)
            sol0=[sol0(1:bcdofs(i)-1); bcvals(i); sol0(bcdofs(i):end)];
        end
        
        toc
        
        fprintf(' p=%d , num_intervals=%d, num_dof = %d\n', p, numix, size(stiff,1))
        disp('Calculating the relative errors in infinity norms...')
        [errw0inf, errw1inf, errw2inf] = calcinfnorm(sol0, elementcounter, p, q, knotu, knotv, b_net, element_nod, coord_ij, nument);
        fprintf('l_inf=%1.15f  w1_inf=%1.15f w2_inf=%1.15f\n', [errw0inf, errw1inf, errw2inf])
        toc
        disp('Calculating the relative errors in square norms...')
        [l2normerr,h1normerr,h2normerr]=calcSqErrNorms(sol0, element_int, elementcounter, p, q, knotu, knotv, b_net, element_nod, coord_ij, nument);
        fprintf('l2=%1.15f  h1=%1.15f h2=%1.15f\n', [l2normerr,h1normerr,h2normerr])
        
        
        toc
    end
end
ploterror(sol0, elementcounter, p, q, knotu, knotv, b_net, lenu, lenv, element_nod, coord_ij)


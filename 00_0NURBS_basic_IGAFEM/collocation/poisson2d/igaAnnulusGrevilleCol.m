% script igaAnnulusGrevilleCol.m
% Solve -Delta(U) + alpha*u = f on a quarter annulus domain with homogeneous Dirichlet boundary
% Greville Collocation method
% INPUT data: p - polynomial order
%             numiset - set that contains the number of elements in each
%             direction


close all
clear all
addpath ../../nurbs-geopdes/inst

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
    
    for numix = numiset % number of intervals along the length
        
        numiy = numix;
        
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
        
        ngaussx = p+1;
        ngaussy = q+1;
        ngaussedge = max(ngaussx,ngaussy);
        
        
        
        %intialize the control net
        coordinates = zeros(numnodes, dim+1);  %allocate one more column for the weights
        coord_ij = zeros(numnodes, dim); %coordinate # -> i,j array
        b_net = zeros(lenu, lenv, dim+1); %same array in the B_{ij} format
        
        %initialize the arrays for the coordinates of the collocations
        %sample points
        
        num_samp = (lenu-2)*(lenv-2);
        sample_points = zeros(num_samp, dim); %coordinates of the sample points in the interior
        %coordinates and orientation of sample points on edges
        sample_pts_left = zeros(lenv, dim+1);
        sample_pts_right = zeros(lenv, dim+1);
        sample_pts_bottom = zeros(lenu, dim+1);
        sample_pts_top = zeros(lenu, dim+1);
        
        sample_counter = 0; %counter for interior sample points
        sample_counter_left = 0;
        sample_counter_right = 0;
        sample_counter_bottom = 0;
        sample_counter_top = 0;
        
        
        samp_span = zeros(num_samp, 2);
        
        for j=1:lenv %for each node in the y direction
            
            %compute the sample points in the x-direction using Greville Abscisae
            coordy = sum(knotv(j+1:j+q))./q;
            
            for i=1:lenu % for each node in the x direction
                index = (j-1)*lenu + i; %the index of the coordinate array
                coordx = sum(knotu(i+1:i+p))./p;
                coordinates(index,:) = [nurbsr.coefs(1,i,j)./nurbsr.coefs(4,i,j), nurbsr.coefs(2,i,j)./nurbsr.coefs(4,i,j), nurbsr.coefs(4,i,j)]; %put the (i,j) node in the coordinate array
                b_net(i,j,:) = coordinates(index,:);
                coord_ij(index,:) = [i, j];
                
                
                if (i>1) && (i<lenu) && (j>1) && (j<lenv)
                    sample_counter = sample_counter +1;
                    sample_points(sample_counter,:) = [coordx, coordy];
                    samp_span(sample_counter,:) = [i,j];
                end
                
                if i==1
                    sample_counter_left = sample_counter_left+1;
                    sample_pts_left(sample_counter_left,:) = [coordx, coordy, 2];
                end
                
                if i==lenu
                    sample_counter_right = sample_counter_right + 1;
                    sample_pts_right(sample_counter_right,:) = [coordx, coordy, 4];
                end
                
                if j==1
                    sample_counter_bottom = sample_counter_bottom + 1;
                    sample_pts_bottom(sample_counter_bottom,:) = [coordx, coordy, 1];
                end
                
                if j==lenv
                    sample_counter_top = sample_counter_top + 1;
                    sample_pts_top(sample_counter_top,:) = [coordx, coordy, 3];
                end
                
            end %for i
        end %for j
        
        luku = length(unique(knotu))-1;
        lukv = length(unique(knotv))-1;
        
        nument = (p+1)*(q+1); %number of nodes for each element
        element_nod = zeros(luku*lukv,nument);
        element_int = zeros(luku*lukv,4);
        
        %loop through each element and compute element_int, element_nod,
        %span_elm
        span_elm = zeros(length(knotu)-1, length(knotv)-1);  %knot-span -> element connectivity
        elementcounter = 0;
        for j=1:length(knotv)-1
            for i=1:length(knotu)-1
                if (abs(knotu(i+1)-knotu(i))>tolcomp) && (abs(knotv(j+1)-knotv(j))>tolcomp)  %the knotspan has non-zero area
                    elementcounter = elementcounter + 1;
                    span_elm(i,j) = elementcounter;
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
        
        %establish sample point -> element connectivity
        samp_elm = zeros(num_samp, 1);
        
        for i=1:num_samp
            curu = sample_points(i,1);
            curv = sample_points(i,2);
            for t1 = 1:length(knotu)-1
                if (curu >= knotu(t1)) && (curu <= knotu(t1+1)) && (knotu(t1) < knotu(t1+1))
                    for t2 = 1:length(knotv)-1
                        if (curv >= knotv(t2)) && (curv <= knotv(t2+1)) && (knotv(t2) < knotv(t2+1))
                            samp_elm(i) = span_elm(t1,t2);
                        end
                    end
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
        
        deriv_order = 2; %we need both derivatives for collocation
        [ M_arr_u ] = makeIndexColBspline( knotu, p, deriv_order);
        [ M_arr_v ] = makeIndexColBspline( knotv, q, deriv_order);
        
        %assembly
        toc
        disp('Assembling linear system...')
        
        %index for storing entries of the LHS
        II = zeros(1, num_samp*nument);
        JJ = zeros(1, num_samp*nument);
        S = zeros(1,  num_samp*nument);
        
        colrhs = zeros(num_samp,1);
        
        for sampcounter=1:num_samp
            
            %retrieve the element and span indexes
            i = samp_elm(sampcounter);
            t1 = samp_span(sampcounter, 1);
            t2 = samp_span(sampcounter, 2);
            
            
            %find the knot span index in each direction
            ni = coord_ij(element_nod(i,end),1);
            nj = coord_ij(element_nod(i,end),2);
            
            %calculate the weights and control points corresponding to the
            %current element
            
            wgts = reshape(b_net(ni-p:ni, nj-q:nj, dim+1), nument, 1);
            cpts = reshape(b_net(ni-p:ni, nj-q:nj, 1:2), nument, dim);
            
            M = M_arr_u(:,:,t1);
            P = M_arr_v(:,:,t2);
            
            
            [R, dR, ddR] = nurbshape2d5(M,P,p,q,wgts);
            
            %calculate the coordinates in the physical space
            coord = R*cpts;
            
            
            %calculate the Jacobian of the transformation
            dxdxi = dR*cpts;
            
            
            % Set up the second derivatives matrix and the matrix of squared first derivatives
            d2xdxi2 = ddR*cpts;
            
            dxdxi2 = [dxdxi(1,1)^2, 2*dxdxi(1,1)*dxdxi(1,2), dxdxi(1,2)^2;...
                dxdxi(1,1)*dxdxi(2,1), dxdxi(1,1)*dxdxi(2,2)+dxdxi(1,2)*dxdxi(2,1), dxdxi(1,2)*dxdxi(2,2);...
                dxdxi(2,1)^2, 2*dxdxi(2,1)*dxdxi(2,2), dxdxi(2,2)^2];
            
            % Solve for first derivatives in global coordinates
            dR = dxdxi\dR;
            
            % Solve for second derivatives in global coordinates
            ddR = dxdxi2\(ddR - d2xdxi2*dR);
            
            loclhsc_delta = ddR(1,:)+ddR(3,:);
            loclhsc_u = R;
            scrt = element_nod(i,:);
            
            [~,~,~,z] = exact_sol(coord(1), coord(2));
            colrhs(sampcounter) = z;
            
            II((sampcounter-1)*nument+1:sampcounter*nument)=sampcounter;
            JJ((sampcounter-1)*nument+1:sampcounter*nument)=scrt;
            S((sampcounter-1)*nument+1:sampcounter*nument)=-loclhsc_delta + alpha*loclhsc_u;
            
        end
        
        collhs = sparse(II,JJ,S,num_samp, numnodes);
        clear II JJ S
        
        toc
        disp('Imposing Dirichlet boundary conditions...')
        [bcdof] = DirichletProj(dirichlet, element_nod, numnodes, knotu, knotv, ngaussedge, nument, b_net, p, q, lenu, lenv, coord_ij, coordinates);
        bcval = zeros(size(bcdof));
        
        
        %impose the Dirichlet boundary conditions
        red_vect = zeros(1, numnodes);
        red_vect(bcdof) = bcval;
        
        colrhs = colrhs - collhs*red_vect';
        collhs(:, bcdof) = [];
        
        toc
        
        disp('Solving the linear system...')
        %condition_number_estimate = condest(collhs)
        sol0 = collhs\colrhs;
        
        [bcdofs, ind] = sort(bcdof);
        bcvals = bcval(ind);
        for i=1:length(bcdof)
            sol0=[sol0(1:bcdofs(i)-1); bcvals(i); sol0(bcdofs(i):end)];
        end
        toc
        fprintf(' p=%d , num_intervals=%d, num_dof = %d\n', p, numix, size(collhs,1))
        
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

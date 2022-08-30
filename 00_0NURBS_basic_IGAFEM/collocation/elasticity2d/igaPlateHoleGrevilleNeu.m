%Script for plate with hole problem using collocation at Greville abscisae
%Neuamnn boundary conditions + symmetry Dirichlet

close all
clear all
addpath ../../nurbs-geopdes/inst  %to use functions from Octave NURBS package


%dimensions of the plate centered at the origin, quarter modelled lies in
%the IInd quadrant.
L = 4;  %length of the plate
rad = 1;  %radius of the hole
tx = 10; %force applied in the x direction
ty = 0;  %force applied in the y direction

[Emod, nu] = materialenu(0.5,0.5);
C=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];

dim = 2; % we are on a two-dimensional plate


for p=5:5  %degree of Nurbs in each direction
    q=p;
    
    for numrefinements = 4:4
        tic
        disp('Pre-processing...')
        
        %material properties
        [E,nu]=materialenu(0.5,0.5);
        
        knotu = [0, 0, 0, 0.5, 1, 1, 1];
        knotv = [0, 0, 0, 1, 1, 1];
        
        weights = [1,(1+1/sqrt(2))/2,(1+1/sqrt(2))/2,1, 1,1,1,1, 1,1,1,1]';
        controlPoints = [-1, 0; -1, sqrt(2)-1; 1-sqrt(2), 1; 0, 1; -2.5, 0; -2.5, 0.75; -0.75, 2.5; 0, 2.5; -4, 0; -4, 4; -4, 4; 0, 4];
        
        %the number of control points in the u and v directions
        lenu = length(knotu)-3;  %number of basis functions in the u direction
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
        
        %number of gauss points for error computation
        ngaussx=p+1;
        ngaussy=q+1;
        
        %the number of control points in the u and v directions
        lenu = length(knotu)-p-1;  %number of basis functions in the u direction
        lenv = length(knotv)-q-1;  %number of basis functions in the v direction
        numnodes = lenu*lenv;      %number of control points
        numix = length(unique(knotu))-1;
        numiy = length(unique(knotv))-1;
        
        
        coordinates = [controlPoints, weights];
        
        num_samp = (lenu-2)*(lenv-2);
        sample_points = zeros(num_samp, dim); %coordinates of the sample points in the interior
        
        coord_ij = zeros(numnodes, dim); %coordinate # -> i,j array
        b_net = zeros(lenu, lenv, dim+1); %same array in the B_{ij} format
        
        
        sample_counter = 0; %counter for interior sample points
        
        %coordinates and orientation of sample points on edges
        sample_pts_left = zeros(lenv, dim+1);
        sample_pts_right = zeros(lenv, dim+1);
        sample_pts_bottom = zeros(lenu, dim+1);
        sample_pts_top = zeros(lenu, dim+1);
        
        sample_counter_left = 0;
        sample_counter_right = 0;
        sample_counter_bottom = 0;
        sample_counter_top = 0;
        
        samp_span = zeros(num_samp, 2);
        for j=1:lenv %for each node in the y direction
            
            %compute the control points in the x-direction using Greville Abscisae
            coordy = sum(knotv(j+1:j+q))./q;
            
            for i=1:lenu % for each node in the x direction
                index = (j-1)*lenu + i; %the index of the coordinate array
                coordx = sum(knotu(i+1:i+p))./p;
                b_net(i,j,:) = coordinates(index,:);
                coord_ij(index,:) = [i, j];
                
                if i==1
                    sample_counter_left = sample_counter_left+1;
                    sample_pts_left(sample_counter_left,:) = [coordx, coordy, 4];
                end
                
                if i==lenu
                    sample_counter_right = sample_counter_right + 1;
                    sample_pts_right(sample_counter_right,:) = [coordx, coordy, 2];
                end
                
                if j==1
                    sample_counter_bottom = sample_counter_bottom + 1;
                    sample_pts_bottom(sample_counter_bottom,:) = [coordx, coordy, 1];
                end
                
                if j==lenv
                    sample_counter_top = sample_counter_top + 1;
                    sample_pts_top(sample_counter_top,:) = [coordx, coordy, 3];
                end
                
                if (i>1) && (i<lenu) && (j>1) && (j<lenv)
                    sample_counter = sample_counter +1;
                    sample_points(sample_counter,:) = [coordx, coordy];
                    samp_span(sample_counter,:) = [i,j];
                end
                
            end %for i
        end %for j
        
        
        
        %initialize the element integration matrix
        %format: [x1][y1][x2][y2]...
        
        luku = length(unique(knotu))-1;
        lukv = length(unique(knotv))-1;
        element_int = zeros(luku*lukv,4);
        nument = (p+1)*(q+1); %number of nodes for each element
        element_nod = zeros(luku*lukv,nument);
        
        %loop through each element and compute element_int, element_nod,
        %element_samp
        span_elm = zeros(length(knotu)-1, length(knotv)-1);  %knot-span -> element connectivity
        elementcounter = 0;
        for j=1:length(knotv)-1
            for i=1:length(knotu)-1
                if (knotu(i+1)>knotu(i) && knotv(j+1)>knotv(j))  %the knotspan has non-zero area
                    %            [i, j]
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
        
        samp_elm_left = zeros(lenv, 1);
        for i=1:lenv
            curu = sample_pts_left(i,1);
            curv = sample_pts_left(i,2);
            for t1 = 1:length(knotu)-1
                if (curu >= knotu(t1)) && (curu <= knotu(t1+1)) && (knotu(t1) < knotu(t1+1))
                    for t2 = 1:length(knotv)-1
                        if (curv >= knotv(t2)) && (curv <= knotv(t2+1)) && (knotv(t2) < knotv(t2+1))
                            samp_elm_left(i) = span_elm(t1,t2);
                        end
                    end
                end
            end
        end
        
        samp_elm_right = zeros(lenv, 1);
        for i=1:lenv
            curu = sample_pts_right(i,1);
            curv = sample_pts_right(i,2);
            for t1 = 1:length(knotu)-1
                if (curu >= knotu(t1)) && (curu <= knotu(t1+1)) && (knotu(t1) < knotu(t1+1))
                    for t2 = 1:length(knotv)-1
                        if (curv >= knotv(t2)) && (curv <= knotv(t2+1)) && (knotv(t2) < knotv(t2+1))
                            samp_elm_right(i) = span_elm(t1,t2);
                        end
                    end
                end
            end
        end
        
        samp_elm_bottom = zeros(lenu, 1);
        for i=1:lenu
            curu = sample_pts_bottom(i,1);
            curv = sample_pts_bottom(i,2);
            for t1 = 1:length(knotu)-1
                if (curu >= knotu(t1)) && (curu <= knotu(t1+1)) && (knotu(t1) < knotu(t1+1))
                    for t2 = 1:length(knotv)-1
                        if (curv >= knotv(t2)) && (curv <= knotv(t2+1)) && (knotv(t2) < knotv(t2+1))
                            samp_elm_bottom(i) = span_elm(t1,t2);
                        end
                    end
                end
            end
        end
        
        samp_elm_top = zeros(lenu, 1);
        for i=1:lenu
            curu = sample_pts_top(i,1);
            curv = sample_pts_top(i,2);
            for t1 = 1:length(knotu)-1
                if (curu >= knotu(t1)) && (curu <= knotu(t1+1)) && (knotu(t1) < knotu(t1+1))
                    for t2 = 1:length(knotv)-1
                        if (curv >= knotv(t2)) && (curv <= knotv(t2+1)) && (knotv(t2) < knotv(t2+1))
                            samp_elm_top(i) = span_elm(t1,t2);
                        end
                    end
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
        
        %         for i=1:length(knotu-1)
        %             line([knotu(i), knotu(i)], [0,1], 'Color', 'black');
        %             hold on
        %         end
        %         for i=1:length(knotv-1)
        %             line([0, 1], [knotv(i), knotv(i)], 'Color', 'black');
        %             hold on
        %         end
        % plot(sample_points(:,1), sample_points(:,2), '.')
        % hold on
        % plot(sample_pts_left(:,1), sample_pts_left(:,2), '.r', sample_pts_right(:,1), sample_pts_right(:,2), '.r')
        % hold on
        % plot(sample_pts_bottom(:,1), sample_pts_bottom(:,2), '.r', sample_pts_top(:,1), sample_pts_top(:,2), '.r')
        %
        %
        
        toc
        disp('Pre-computing B-Splines...')
        
        %loop through the knot spans in each direction and pre-compute
        %Bsplines and connectivity arrays
        
        deriv_order = 2; %we need both derivatives for collocation
        [ M_arr_u ] = makeIndexColBspline( knotu, p, deriv_order);
        [ M_arr_v ] = makeIndexColBspline( knotv, q, deriv_order);
        
        toc
        %assembly
        disp('Assembling linear system...')
        %  collhs = sparse(dim*num_samp, dim*numnodes);
        colrhs = zeros(dim*num_samp,1);
        
        
        %index for storing entries of the LHS
        II = zeros(1, 2*dim*num_samp*nument);
        JJ = zeros(1, 2*dim*num_samp*nument);
        S = zeros(1,  2*dim*num_samp*nument);
        
        
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
            
            
            loclhsc = zeros(2, 2*(1+p)*(1+q));
            scrt = zeros(1, 2*(1+p)*(1+q));
            locrhsc = zeros(2,1);
            
            %vectorized global matrix assembly
            mul_constant = Emod/(1-nu^2);
            loclhsc(1, 1:2:end-1) = ddR(1,:) + (1-nu)/2*ddR(3,:);
            loclhsc(2, 1:2:end-1) = (1+nu)/2*ddR(2,:);
            loclhsc(1, 2:2:end) = (1+nu)/2*ddR(2,:);
            loclhsc(2, 2:2:end) = ddR(3,:)+(1-nu)/2*ddR(1,:);
            loclhsc = mul_constant*loclhsc;
            scrt(1:2:end-1) = 2*element_nod(i,:)-1;
            scrt(2:2:end) = 2*element_nod(i,:);
            
            %evaluate right-hand side
            [z1, z2] = f(coord(1), coord(2));
            colrhs(2*sampcounter-1:2*sampcounter) = [z1;z2];
            
            
            II(4*(sampcounter-1)*nument+1:(4*sampcounter-2)*nument)=2*sampcounter-1;
            II((4*sampcounter-2)*nument+1: 4*sampcounter*nument) = 2*sampcounter;
            JJ(4*(sampcounter-1)*nument+1: 4*sampcounter*nument) = [scrt, scrt];
            S(4*(sampcounter-1)*nument+1:(4*sampcounter-2)*nument) = loclhsc(1,:);
            S((4*sampcounter-2)*nument+1: 4*sampcounter*nument) = loclhsc(2,:);
            
        end
        
        collhs = sparse(II, JJ, S, dim*num_samp, dim*numnodes);
        clear II JJ S
        
        
        %Neumann boundary
        toc
        disp('Imposing Neumann boundary conditions...')
        
        neumann_samp_points = [sample_pts_bottom(2:end-1,:); sample_pts_right; sample_pts_top(end-1:-1:2,:); sample_pts_left(end:-1:1,:)];
        
        %plot the neumann points
        plot(neumann_samp_points(:,1), neumann_samp_points(:,2), '.r')
        hold on
        
        neumann_samp_elm = [samp_elm_bottom(2:end-1,:); samp_elm_right; samp_elm_top(end-1:-1:2,:); samp_elm_left(end:-1:1,:)];
        num_neumann_samp = size(neumann_samp_points, 1);
        
        bot_samp = [2*(lenu+lenv)-4, 1:lenu-1];
        top_samp = lenu-2+lenv:2*lenu-3+lenv;
        left_samp = 2*lenu-3+lenv:2*(lenu+lenv)-4;
        right_samp = lenu-1:lenu+lenv-2;
        
        
        
        
        % neumann_samp_points = zeros(num_neumann_samp, dim);
        neulhs = sparse(2*num_neumann_samp, dim*numnodes);
        neurhs = zeros(2*num_neumann_samp, 1);
        
        
        %loop over the elements
        for i=1:size(neumann_samp_points, 1)
            
            orient = neumann_samp_points(i,3);
            %calculate the element index
            elm_index = neumann_samp_elm(i);
            ucoord = neumann_samp_points(i,1);
            vcoord = neumann_samp_points(i,2);
            
            [~, dN, ~, coord, normal] = nurbshaped(elm_index, ucoord, vcoord, knotu, knotv, b_net, p, q, lenu, lenv,...
                element_nod, coord_ij, orient);
            
            l =normal(1);
            m =normal(2);
            
            loclhsc = zeros(2, 2*(1+p)*(1+q));
            scrt = zeros(1, 2*(1+p)*(1+q));
            locrhsc = zeros(2,1);
            
            for t=1:(1+p)*(1+q)
                
                loclhsc(1:2, 2*t-1:2*t) = Emod/(1-nu^2).*[l*dN(1,t) + m*(1-nu)/2*dN(2,t),...
                    l*nu*dN(2,t)+m*(1-nu)/2*dN(1,t); m*nu*dN(1,t)+l*(1-nu)/2*dN(2,t),...
                    m*dN(2,t)+l*(1-nu)/2*dN(1,t)];
                scrt(2*t-1:2*t) = [2*element_nod(elm_index,t)-1, 2*element_nod(elm_index,t)];
                
            end
            
            [stress] = ghole(coord(1), coord(2), rad, tx);
            
            taux = normal(1)*stress(1) + normal(2)*stress(3);
            tauy = normal(1)*stress(3) + normal(2)*stress(2);
            
            neurhs(2*i-1:2*i) = [taux;tauy];
            neulhs(2*i-1:2*i,scrt) = loclhsc;
        end
        
        toc
        disp('Imposing Dirichlet boundary conditions...')
        plot(sample_points(:,1), sample_points(:,2), '.')
        hold on
        
        [dirichlet_nodesl] = DirichletProjCo(leftedge, element_nod, p, q);
        [dirichlet_nodesr] = DirichletProjCo(rightedge, element_nod, p, q);
        [dirichlet_nodesb] = DirichletProjCo(bottomedge, element_nod, p, q);
        [dirichlet_nodest] = DirichletProjCo(topedge, element_nod, p, q);
        
        num_left_nodes = length(dirichlet_nodesl);
        num_right_nodes = length(dirichlet_nodesr);
        num_bottom_nodes = length(dirichlet_nodesb);
        num_top_nodes = length(dirichlet_nodest);
        
        bcsaml = zeros(1*(num_left_nodes),1);
        bcdofl = zeros((num_left_nodes),1);
        for i=1:num_left_nodes  %for each node on the left edge
            %calculate global node index for y displacement
            bcdofl(i) = 2*dirichlet_nodesl(i);
            bcsaml(i) = 2*left_samp(end+1-i);
        end
        
        bcdofr = zeros((num_right_nodes),1);
        bcsamr = zeros(1*(num_right_nodes),1);
        for i=1:num_right_nodes %for each node on the right edge
            %calculate global node index
            bcdofr(i) = 2*dirichlet_nodesr(i)-1;
            bcsamr(i) = 2*right_samp(i)-1;
        end
        
        
        bcdof = [bcdofl; bcdofr];
        bcval = zeros(size(bcdof));
        bcsam = [bcsaml; bcsamr];
        
        
%        [ndlhs, ndrhs] = feaplycloc2(neulhs, neurhs, bcdof, bcsam, bcval);      
         ndlhs = neulhs;
         ndrhs = neurhs;
         ndlhs(bcsam,:) = [];
         ndrhs(bcsam,:) = [];
%         size(ndlhs)
        collhs = [collhs; ndlhs];
        colrhs = [colrhs; ndrhs];
        
        
        %impose the Dirichlet boundary conditions
        red_vect = zeros(1, dim*numnodes);
        red_vect(bcdof) = bcval;
        
        
        colrhs = colrhs - collhs*red_vect';
        collhs(:, bcdof) = [];
        
        toc
        disp('Solving the linear system...')
        size(collhs)
        x=collhs\colrhs;
        
        [bcdofs, ind] = sort(bcdof);
        bcvals = bcval(ind);
        for i=1:length(bcdof)
            x=[x(1:bcdofs(i)-1); bcvals(i); x(bcdofs(i):end)];
        end
        toc
        
        disp('Computing the error...')
        
        [l2errrel, relenergynorm] = h1err(x, element_nod, ngaussx, ngaussy, elementcounter, element_int, coord_ij, knotu, knotv, b_net, p, q, C, rad, Emod, nu, tx);
        fprintf('p= %d num_ref = %d dof = %d l2err = %1.15f energy_err= %1.15f \n', p, numrefinements, length(x), l2errrel, relenergynorm)
        
        toc
    end
end

disp('Visualizing the result...')
visualize
toc

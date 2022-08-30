% Solve -Delta(U) + alpha*u = f on hollow sphere domain with homogeneous Dirichlet boundary
% Greville method 
close all
clear all

addpath ../../nurbs-geopdes/inst
numiset = 10:5:10; % number of intervals set
[~, ~, ~, ~, ~, alpha] = exact_sol(0,0,0); %solve -Delta u + alpha*y = f.

indexp = 0;
for p=3:3 %polynomial degree
     
    q = p;
    r = p;
    indexp = indexp+1;
    indexElement = 0;
    for numix = numiset % number of intervals along the length
        
        indexElement = indexElement+1;
        numiy = numix;
        numiz = numix;
        
        tic
        disp(datestr(now))
        disp('Pre-processing...')
        
        knotu = [0, 0, 0, 1, 1, 1];
        knotv = [0, 0, 0, 1, 1, 1];
        knotw = [0, 0, 1, 1];
        %knotw = [0, 0, 0, 1, 1, 1];
        
        h=1;
        
        
        weights = [1, 1/sqrt(2), 1, 1/sqrt(2), 1/2, 1/sqrt(2), 1, 1/sqrt(2), 1, 1, 1/sqrt(2), 1, 1/sqrt(2), 1/2, 1/sqrt(2), 1, 1/sqrt(2), 1]';
        controlPoints = [10-h,0,0; 10-h, 0, 10-h; 0, 0, 10-h; 10-h, 10-h, 0; 10-h, 10-h, 10-h; 0, 0, 10-h; 0, 10-h, 0; 0, 10-h, 10-h; 0, 0, 10-h; ...
            10+h,0,0; 10+h, 0, 10+h; 0, 0, 10+h; 10+h, 10+h, 0; 10+h, 10+h, 10+h; 0, 0, 10+h; 0, 10+h, 0; 0, 10+h, 10+h; 0, 0, 10+h];
        
        %the number of control points in the u and v directions
        lenu = 3;  %number of basis functions in the u direction
        lenv = 3;  %number of basis functions in the v direction
        lenw = 2;  %number of basis functions in the w direction
        coordinates = [controlPoints, weights];
        
        %convert control points/weights to coefs format in NURBS toolbox
        coefs = zeros(4,lenu,lenv);
        index = 0;
        for k=1:lenw
            for i=1:lenu
                for j=1:lenv
                    index = index + 1;
                    xcoord = coordinates(index,1);
                    ycoord = coordinates(index,2);
                    zcoord = coordinates(index,3);
                    wght = coordinates(index,4);
                    coefs(1,i,j,k) = xcoord*wght;
                    coefs(2,i,j,k) = ycoord*wght;
                    coefs(3,i,j,k) = zcoord*wght;
                    coefs(4,i,j,k) = wght;
                end
            end
        end
        knots =  {knotu, knotv, knotw};
        srf = nrbmak(coefs, knots); %generates the coarsest nurbs object, p=2, q=2, r=1
        srf = nrbdegelev(srf, [p,q,r]-(srf.order-1)); %elevate to given p,q,r
        
        %refine by knot insertion
        knot1 = linspace(0,1,numix+1);
        knot2 = linspace(0,1,numiy+1);
        knot3 = linspace(0,1,numiz+1);
        %knot1 =[];
        %knot2 = [];
        srf = nrbkntins(srf, {knot1(2:end-1) knot2(2:end-1) knot3(2:end-1)});
              
        %convert NURBS toolbox format back to control points/weights
        knotu = srf.knots{1};
        knotv = srf.knots{2};
        knotw = srf.knots{3};
        
        %the number of control points in the u and v directions
        lenu = length(knotu)-p-1;  %number of basis functions in the u direction
        lenv = length(knotv)-q-1;  %number of basis functions in the v direction
        lenw = length(knotw)-r-1;  %number of basis functions in the w direction
        
        
        dim = 3; % we are on a three-dimensional world
        b_net = zeros(lenu, lenv, lenw, dim+1);
        coord_ijk = zeros(lenu*lenv*lenw, dim);
        coordinates = zeros(lenu*lenv*lenw, dim+1);
        index = 0;
        for k=1:lenw
            for j=1:lenv
                for i=1:lenu
                    index = index + 1;
                    xcoord = srf.coefs(1,i,j,k)/srf.coefs(4,i,j,k);
                    ycoord = srf.coefs(2,i,j,k)/srf.coefs(4,i,j,k);
                    zcoord = srf.coefs(3,i,j,k)/srf.coefs(4,i,j,k);
                    wght = srf.coefs(4,i,j,k);
                    b_net(i,j,k,:) = [xcoord, ycoord, zcoord, wght];
                    coord_ijk(index,:) = [i, j, k];
                    coordinates(index, :) = [xcoord, ycoord, zcoord, wght];
                    
                end
            end
        end
        numnodes = lenu*lenv*lenw;      %number of control points
        
        ngaussx = p+1;
        ngaussy = q+1;
        ngaussz = r+1;
        ngaussedge = max([ngaussx,ngaussy,ngaussz]);
        
        %initialize the arrays for the coordinates of the collocations
        %sample points
        
        num_samp = (lenu-2)*(lenv-2)*(lenw-2);
        sample_points = zeros(num_samp, dim); %coordinates of the sample points in the interior
        
        sample_counter = 0; %counter for interior sample points
        samp_span = zeros(num_samp, dim);
        
        
        for k=1:lenw
            %compute the control points in the z-direction using Greville Abscisae
            coordz = sum(knotw(k+1:k+r))./r;
            for j=1:lenv %for each node in the y direction
                %compute the control points in the y-direction using Greville Abscisae
                coordy = sum(knotv(j+1:j+q))./q;
                for i=1:lenu % for each node in the x direction
                    
                    index = index + 1;
                    %compute the control points in the x-direction using Greville Abscisae
                    coordx = sum(knotu(i+1:i+p))./p;
                    
                    % coordinates(index,:) = [coordx, coordy, coordz, 1]; %put the (i,j,k) node in the coordinate array with weight 1
                    % b_net(i,j,k,:) = coordinates(index,:);
                    % coord_ijk(index,:) = [i, j, k];
                    if (i>1) && (i<lenu) && (j>1) && (j<lenu) && (k>1) && (k<lenu) %if we are at an interior control point
                        sample_counter = sample_counter +1;  %set it as sample point for collocation
                        sample_points(sample_counter,:) = [coordx, coordy, coordz];
                        samp_span(sample_counter,:) = [i,j,k];
                    end
                end %for i
            end %for j
        end %for k
        
        
        luku = length(unique(knotu))-1;
        lukv = length(unique(knotv))-1;
        lukw = length(unique(knotw))-1;
        
        nument = (p+1)*(q+1)*(r+1); %number of nodes for each element
        element_nod = zeros(luku*lukv*lukw,nument);
        element_int = zeros(luku*lukv*lukw,6);
        
        
        %loop through each element and compute element_int, element_nod,
        %span_elm
        span_elm = zeros(length(knotu)-1, length(knotv)-1, length(knotw)-1);  %knot-span -> element connectivity
        elementcounter = 0;
        for k=1:length(knotw)-1
            for j=1:length(knotv)-1
                for i=1:length(knotu)-1
                    if (knotu(i+1)>knotu(i)) && (knotv(j+1) >knotv(j)) && (knotw(k+1) > knotw(k))  %the knotspan has non-zero area
                        elementcounter = elementcounter + 1;
                        span_elm(i,j,k) = elementcounter;
                        element_int(elementcounter,:) = [knotu(i), knotv(j), knotw(k), knotu(i+1), knotv(j+1),knotw(k+1)];
                        tcount = 0;
                        currow = zeros(1, nument);
                        %now we add the nodes from i-p...i in the u
                        %direction, j-q...j in the v direction, k-r...k in
                        %w direction
                        for t3 = k-r:k
                            for t2=j-q:j
                                for t1 = i-p:i
                                    tcount = tcount + 1;
                                    currow(tcount) = t1+(t2-1)*lenu+(t3-1)*lenu*lenv;
                                end
                            end
                        end
                        element_nod(elementcounter,:)=currow;
                    end
                end
            end
        end
        
        %establish sample point -> element connectivity
        samp_elm = zeros(num_samp, 1);
        for i=1:num_samp
            curu = sample_points(i,1);
            curv = sample_points(i,2);
            curw = sample_points(i,3);
            for t1 = 1:length(knotu)-1
                if (curu >= knotu(t1)) && (curu <= knotu(t1+1)) && (knotu(t1) < knotu(t1+1))
                    for t2 = 1:length(knotv)-1
                        if (curv >= knotv(t2)) && (curv <= knotv(t2+1)) && (knotv(t2) < knotv(t2+1))
                            for t3 = 1:length(knotw)-1
                                if (curw >= knotw(t3)) && (curw <= knotw(t3+1)) && (knotw(t3) < knotw(t3+1))
                                    samp_elm(i) = span_elm(t1,t2,t3);
                                end
                            end
                        end
                    end
                end
            end
        end
               
        
        %identify the nodes on each face to impose Dirichlet boundary cond.
        top_nodes = lenu*lenv*(lenw-1)+1:lenu*lenv*lenw;
        bottom_nodes = 1:lenu*lenv;
        front_nodes = zeros(1, lenw*lenu);
        back_nodes = zeros(1, lenw*lenu);
        left_nodes = zeros(1, lenw*lenv);
        right_nodes = zeros(1, lenw*lenv);
        for i=1:lenw
            front_nodes((i-1)*lenu+1:i*lenu) = (lenu*lenv)*(i-1)+1:(lenu*lenv)*(i-1)+lenu;
            back_nodes((i-1)*lenu+1:i*lenu) = i*(lenu*lenv):-1:i*(lenu*lenv)-lenu+1;
            left_nodes((i-1)*lenv+1:i*lenv) = (i-1)*lenu*lenv+1:lenu:(i-1)*lenu*lenv+lenu*(lenv-1)+1;
            right_nodes((i-1)*lenv+1:i*lenv) = (i-1)*lenu*lenv+lenu:lenu:i*lenu*lenv;
        end
        
        %pick which nodes correspond to neumann and dirichlet conditions
        
        neumann = [];
        dirichlet = unique([top_nodes, bottom_nodes, front_nodes, back_nodes, left_nodes, right_nodes]);
        
        toc
        disp('Pre-computing B-Splines...')
        
        %loop through the knot spans in each direction and pre-compute
        %Bsplines and connectivity arrays
        
        deriv_order = 2; %we need 2nd derivatives for collocation
        [ M_arr_u ] = makeIndexColBspline( knotu, p, deriv_order);
        [ M_arr_v ] = makeIndexColBspline( knotv, q, deriv_order);
        [ M_arr_w ] = makeIndexColBspline( knotw, r, deriv_order);
        
        
        %assembly
        toc
        disp('Assembling linear system...')
        
        %index for storing entries of the LHS
        II = zeros(1, num_samp*nument);
        JJ = zeros(1, num_samp*nument);
        S = zeros(1, num_samp*nument);
        colrhs = zeros(num_samp,1);
        
        for sampcounter=1:num_samp
            i = samp_elm(sampcounter);
            t1 = samp_span(sampcounter, 1);
            t2 = samp_span(sampcounter, 2);
            t3 = samp_span(sampcounter, 3);
            
            
            %find the knot span index in each direction
            ni = coord_ijk(element_nod(i,end),1);
            nj = coord_ijk(element_nod(i,end),2);
            nk = coord_ijk(element_nod(i,end),3);
            
            %calculate the weights and control points corresponding to the
            %current element
            
            wgts = reshape(b_net(ni-p:ni, nj-q:nj, nk-r:nk, dim+1), nument, 1);
            cpts = reshape(b_net(ni-p:ni, nj-q:nj, nk-r:nk, 1:3), nument, 3);
            
            M = M_arr_u(:,:,t1);
            P = M_arr_v(:,:,t2);
            Q = M_arr_w(:,:,t3);
            
            [R, dR, ddR] = nurbshape3d5(M,P,Q,p,q,r,wgts);
            
            % calculate coordinates in physical space
            coord = R*cpts;
            
            % Compute Jacobian matrix
            dxdxi = dR*cpts;           
                        
            % Set up the second derivatives matrix and the matrix of squared first derivatives
            d2xdxi2 = ddR*cpts;
            
            dxdxi2 = [ dxdxi(1,1)^2   dxdxi(1,2)^2  dxdxi(1,3)^2    2*dxdxi(1,1)*dxdxi(1,2)  2*dxdxi(1,1)*dxdxi(1,3)  2*dxdxi(1,2)*dxdxi(1,3);
                dxdxi(2,1)^2   dxdxi(2,2)^2  dxdxi(2,3)^2    2*dxdxi(2,1)*dxdxi(2,2)  2*dxdxi(2,1)*dxdxi(2,3)  2*dxdxi(2,2)*dxdxi(2,3);
                dxdxi(3,1)^2   dxdxi(3,2)^2  dxdxi(3,3)^2    2*dxdxi(3,1)*dxdxi(3,2)  2*dxdxi(3,1)*dxdxi(3,3)  2*dxdxi(3,2)*dxdxi(3,3);
                dxdxi(1,1)*dxdxi(2,1)   dxdxi(1,2)*dxdxi(2,2)   dxdxi(1,3)*dxdxi(2,3)  dxdxi(1,1)*dxdxi(2,2)+dxdxi(2,1)*dxdxi(1,2)   dxdxi(1,1)*dxdxi(2,3)+dxdxi(2,1)*dxdxi(1,3)   dxdxi(1,2)*dxdxi(2,3)+dxdxi(2,2)*dxdxi(1,3);
                dxdxi(1,1)*dxdxi(3,1)   dxdxi(1,2)*dxdxi(3,2)   dxdxi(1,3)*dxdxi(3,3)  dxdxi(1,1)*dxdxi(3,2)+dxdxi(3,1)*dxdxi(1,2)   dxdxi(1,1)*dxdxi(3,3)+dxdxi(3,1)*dxdxi(1,3)   dxdxi(1,2)*dxdxi(3,3)+dxdxi(3,2)*dxdxi(1,3);
                dxdxi(2,1)*dxdxi(3,1)   dxdxi(2,2)*dxdxi(3,2)   dxdxi(2,3)*dxdxi(3,3)  dxdxi(2,1)*dxdxi(3,2)+dxdxi(3,1)*dxdxi(2,2)   dxdxi(2,1)*dxdxi(3,3)+dxdxi(3,1)*dxdxi(2,3)   dxdxi(2,2)*dxdxi(3,3)+dxdxi(3,2)*dxdxi(2,3) ];
            
            % Solve for first derivatives in global coordinates
            dR = dxdxi\dR;
            
            % Solve for second derivatives in global coordinates
            ddR = dxdxi2\(ddR - d2xdxi2*dR);
            
            loclhsc_delta = ddR(1,:)+ddR(2,:)+ddR(3,:);
            loclhsc_u = R;
            scrt = element_nod(i,:);
            
            [~, ~, ~, ~, f] = exact_sol(coord(1),coord(2),coord(3));
            II((sampcounter-1)*nument+1:sampcounter*nument)=sampcounter;
            JJ((sampcounter-1)*nument+1:sampcounter*nument)=scrt;
            S((sampcounter-1)*nument+1:sampcounter*nument)=-loclhsc_delta + alpha*loclhsc_u;
            colrhs(sampcounter) = f;
            
        end
        collhs = sparse(II,JJ,S,num_samp, numnodes);
        clear II JJ S
        
        toc
        disp('Imposing Dirichlet boundary conditions...')
        bcdof = dirichlet;
        bcval = zeros(size(bcdof));
        
        red_vect = zeros(1, numnodes);
        red_vect(bcdof) = bcval;
        
        colrhs = colrhs - collhs*red_vect';
        collhs(:, bcdof) = [];
        
        toc
        
        disp('Preconditioning')
        [L_pre, U_pre] = ilu(collhs);
         
        toc
        disp('Solving the linear system...')
       
        [sol0,fl2,rr2,it2,rv2] = gmres(collhs,colrhs,[],1e-13,1000,L_pre,U_pre);
        fprintf('GMRES exited with flag %d\n', fl2)
        fprintf('Residual value: %1.15g\n', rr2)
        fprintf('Number of iterations: %d\n', it2(2))
        
        disp('Post-processing... ')
        
        [bcdofs, ind] = sort(bcdof);
        bcvals = bcval(ind);
        for i=1:length(bcdof)
            sol0=[sol0(1:bcdofs(i)-1); bcvals(i); sol0(bcdofs(i):end)];
        end
        toc
        fprintf(' p=%d , num_intervals=%d, num_dof = %d\n', p, numix, size(collhs,1))
            
        disp('Calculating the relative errors in square norms...')        
        [l2normerr,h1normerr,h2normerr]=calcSqErrNorms(sol0, element_int, elementcounter, p, q, r, knotu, knotv, knotw, b_net, element_nod, coord_ijk);
        sprintf('%1.15f  %1.15f %1.15f\n', [l2normerr,h1normerr,h2normerr])
              
        toc
    end
end
%
disp('Plotting the error...')
ploterr(sol0, span_elm, p, q, r, knotu, knotv, knotw, b_net, element_nod, coord_ijk,lenu,lenv,lenw);

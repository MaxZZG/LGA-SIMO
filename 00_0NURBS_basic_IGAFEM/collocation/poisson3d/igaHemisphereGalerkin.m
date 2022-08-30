% Solve -Delta(U) + alpha*u = f on a hollow sphere with homogeneous Dirichlet boundary
% Galerkin method

close all
clear all

addpath ../../nurbs-geopdes/inst

numiset = 10:5:10; % number of intervals set
[~, ~, ~, ~, ~, alpha] = exact_sol(0,0,0); %solve -Delta u + alpha*y = f.

indexp = 0;
for p=3:3  %polynomial degree

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
        
        %figure
        %nrbkntplot(srf)
        %nrbctrlplot(srf)
        
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
                    
                end
            end
        end
        numnodes = lenu*lenv*lenw;      %number of control points
        
        ngaussx = p+1;
        ngaussy = q+1;
        ngaussz = r+1;
        ngaussedge = max([ngaussx,ngaussy,ngaussz]);
        
        luku = length(unique(knotu))-1;
        lukv = length(unique(knotv))-1;
        lukw = length(unique(knotw))-1;
        
        nument = (p+1)*(q+1)*(r+1); %number of nodes for each element
        element_nod = zeros(luku*lukv*lukw,nument);
        element_int = zeros(luku*lukv*lukw,6);
        
        %loop through each element and compute element_int, element_nod
        span_elm = zeros(length(knotu)-1, length(knotv)-1, length(knotw)-1);  %knot-span -> element connectivity
        elementcounter = 0;
        for k=1:length(knotw)-1
            for j=1:length(knotv)-1
                for i=1:length(knotu)-1
                    if (knotu(i+1)>knotu(i)) && (knotv(j+1) >knotv(j)) && (knotw(k+1) > knotw(k))  %the knotspan has non-zero area
                        %            [i, j]
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
        
        
        
        %pick wich edges correspond to neumann and dirichlet conditions
        
        neumann = [];
        dirichlet = unique([top_nodes, bottom_nodes, front_nodes, back_nodes, left_nodes, right_nodes]);
        
        
        toc
        disp('Pre-computing B-Splines...')
        
        %loop through the knot spans in each direction and pre-compute
        %Bsplines and connectivity arrays
        
        [gpx,gwx]=genGP_GW(ngaussx);
        [gpy,gwy]=genGP_GW(ngaussy);
        [gpz,gwz]=genGP_GW(ngaussz);
        
        deriv_order = 1; %we only need 1st derivatives for Galerkin method
        [ pt_index_u, M_arr_u ] = makeIndexBspline( knotu, p, deriv_order, gpx );
        [ pt_index_v, M_arr_v ] = makeIndexBspline( knotv, q, deriv_order, gpy );
        [ pt_index_w, M_arr_w ] = makeIndexBspline( knotw, r, deriv_order, gpz );
        
        %assembly
        toc
        disp('Assembling linear system...')
        element_block_size = numix;
        II = zeros(1, element_block_size*nument^2);
        JJ = zeros(1, element_block_size*nument^2);
        indexcounter = 0;
        
        S = zeros(1, element_block_size*nument^2);
        
        stiff = sparse(numnodes, numnodes);
        %mass = sparse(numnodes, numnodes);
        rhs = zeros(numnodes, 1);
        
        for i=1:elementcounter
            
            %for each element calculate the value of the basis
            %functions
            
            ximin = element_int(i,1);
            etamin = element_int(i,2);
            zetamin = element_int(i,3);
            ximax = element_int(i,4);
            etamax = element_int(i,5);
            zetamax = element_int(i,6);
            
            scalefac = (ximax - ximin)*(etamax - etamin)*(zetamax-zetamin)/8;
            
            
            localstiff = zeros(nument, nument); %local stiffness for laplacian term
            localmass = zeros(nument, nument); %local mass matrix
            localrhs = zeros(nument, 1);
            scrtx = element_nod(i,:);
            
            %find the knot span index in each direction
            ni = coord_ijk(element_nod(i,end),1);
            nj = coord_ijk(element_nod(i,end),2);
            nk = coord_ijk(element_nod(i,end),3);
            
            %calculate the weights and control points corresponding to the
            %current element
            
            wgts = reshape(b_net(ni-p:ni, nj-q:nj, nk-r:nk, dim+1), nument, 1);
            cpts = reshape(b_net(ni-p:ni, nj-q:nj, nk-r:nk, 1:3), nument, 3);
            
            %for each shape function in shapelist compute the matrices B, E
            for ii=1:ngaussx
                for jj=1:ngaussy
                    for kk=1:ngaussz
                        
                        glob_index_u = pt_index_u(ni, ii);
                        glob_index_v = pt_index_v(nj, jj);
                        glob_index_w = pt_index_w(nk, kk);
                        
                        M = M_arr_u(:,:,glob_index_u);
                        P = M_arr_v(:,:,glob_index_v);
                        Q = M_arr_w(:,:,glob_index_w);
                        
                        [R, dRdx] = nurbshape3d4(M,P,Q,p,q,r,wgts);
                        
                        %calculate the coordinates in the physical space
                        coord = R*cpts;
                        
                        %calculate the Jacobian of the transformation
                        dxdxi = dRdx*cpts;
                        J = det(dxdxi);
                        
                        %calculate the derivatives with respect to physical
                        %space
                        dRdx = dxdxi\dRdx;
                        
                        [~, ~, ~, ~, f] = exact_sol(coord(1),coord(2),coord(3));
                        tic
                        localrhs = localrhs + f.*R'.*J.*scalefac.*gwx(ii).*gwy(jj).*gwz(kk);
                        if alpha == 0
                            localstiff = localstiff + dRdx'*dRdx.*J.*scalefac.*gwx(ii).*gwy(jj).*gwz(kk);
                        else
                            localstiff = localstiff + (dRdx'*dRdx+alpha*(R'*R)).*J.*scalefac.*gwx(ii).*gwy(jj).*gwz(kk);
                        end
                        
                    end
                end
            end
            if (indexcounter+nument^2)>element_block_size*nument^2
                stiff = stiff + sparse(II,JJ,S,numnodes,numnodes);
                indexcounter = 0;                
            end
            II(indexcounter+1:indexcounter+nument^2) = repmat(scrtx,1,nument);
            JJ(indexcounter+1:indexcounter+nument^2) = reshape(repmat(scrtx',1,nument)',1,nument^2);
            S(indexcounter+1:indexcounter+nument^2) = reshape(localstiff,1,nument^2);
            indexcounter = indexcounter + nument^2;
            
            
            rhs(scrtx) = rhs(scrtx) + localrhs;
        end
        
        stiff = stiff + sparse(II,JJ,S,numnodes,numnodes);
        
        clear II JJ S
        
        
        toc
        disp('Imposing Dirichlet boundary conditions...')
        [bcdof] = dirichlet;
        bcval = zeros(size(bcdof));
        [stiff,rhs] = feaplyc2_poisson(stiff,rhs,bcdof,bcval);
        
        toc
        disp('Preconditioning...')
      
        L_pre = ichol(stiff);
        
        
        toc
        disp('Solving the linear system...')
        [sol0,fl2,rr2,it2,rv2] = pcg(stiff,rhs,1e-14,200,L_pre,L_pre');
        fprintf('PCG exited with flag %d\n', fl2)
        fprintf('Residual value: %1.15g\n', rr2)
        fprintf('Number of iterations: %d\n', it2)
        szstiff = size(stiff);
        clear stiff L_pre
        
        
        toc
        
        disp('Post-processing... ')
        
        [bcdofs, ind] = sort(bcdof);
        bcvals = bcval(ind);
        for i=1:length(bcdof)
            sol0=[sol0(1:bcdofs(i)-1); bcvals(i); sol0(bcdofs(i):end)];
        end
        
        
        fprintf(' p=%d , num_intervals=%d, num_dof = %d\n', p, numix, szstiff(1))
        
             
        disp('Calculating the relative errors in square norms...')
        [l2normerr,h1normerr,h2normerr]=calcSqErrNorms(sol0, element_int, elementcounter, p, q, r, knotu, knotv, knotw, b_net, element_nod, coord_ijk);
        fprintf('%1.15f  %1.15f %1.15f\n', [l2normerr,h1normerr,h2normerr])
       
        toc
    end
end

disp('Plotting the error on the midsurface...')
ploterr(sol0, span_elm, p, q, r, knotu, knotv, knotw, b_net, element_nod, coord_ijk,lenu,lenv,lenw);

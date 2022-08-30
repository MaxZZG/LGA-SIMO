function [l2errrel, relenergynorm] = h1err(x, element_nod, ngaussx, ngaussy, elementcounter, element_int, coord_ij, knotu, knotv, b_net, p, q, C, rad, Emod, nu, tx)

l2errvar = 0; %l2 norm of the error
l2norm = 0; %l2 norm of the solution
energynorm = 0; %energy norm of the computed solution
energyerrnorm = 0; %energy norm of the error
dim = 2;
nument = (p+1)*(q+1);
invC = inv(C);


%loop through the knot spans in each direction and pre-compute
%Bsplines and connectivity arrays

[gpx,gwx]=genGP_GW(ngaussx);
[gpy,gwy]=genGP_GW(ngaussy);

deriv_order = 1; %we only need 1st derivatives for computing the errors
[ pt_index_u, M_arr_u ] = makeIndexBspline( knotu, p, deriv_order, gpx );
[ pt_index_v, M_arr_v ] = makeIndexBspline( knotv, q, deriv_order, gpy );


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
    
    
    scrtx = element_nod(i,:);
    
    %calculate the weights and control points corresponding to the
    %current element
    
    wgts = reshape(b_net(ni-p:ni, nj-q:nj, dim+1), nument, 1);
    cpts = reshape(b_net(ni-p:ni, nj-q:nj, 1:2), nument, dim);
    
    
    dispmatx = zeros(ngaussy,ngaussx);
    dispmaty = zeros(ngaussy,ngaussx);
    exdispmatx = zeros(ngaussy, ngaussx);
    exdispmaty = zeros(ngaussy, ngaussx);
    
    temp = 0; %temporary variable for norm computation
    
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
            
            %evaluate the exact displacements
            eldisp = holeu_d([coord(1), coord(2)], rad, Emod, nu, tx);
            exdispmatx(jj, ii) = eldisp(1);
            exdispmaty(jj, ii) = eldisp(2);
            
            %evaluate the exact stresses
            estress0 = ghole(coord(1), coord(2), rad, tx);
            
            %calculate the Jacobian of the transformation
            dxdxi = dRdx*cpts;
            J = det(dxdxi);
            
            %calculate the derivatives with respect to physical
            %space
            dRdx = dxdxi\dRdx;
            dRdx = dRdx';
            
            %double the entries in scrtx
            dscrtx = reshape([2*scrtx-1; 2*scrtx],1,2*nument);
            
            B(1:2:2*nument-1,1) = dRdx(:,1);
            B(2:2:2*nument,2) = dRdx(:,2);
            B(1:2:2*nument-1,3) = dRdx(:,2);
            B(2:2:2*nument,3) = dRdx(:,1);
            estress = C*B'*x(dscrtx);
            
            %evaluate the computed displacements
            dispmatx(jj,ii) = dispmatx(jj,ii) + R*x(2*scrtx-1);
            dispmaty(jj,ii) = dispmaty(jj,ii) + R*x(2*scrtx);
            
            
            %integrate to evaluate the norms
            temp = temp + (estress0-estress)'*invC*(estress0-estress)*J.*scalefac.*gwx(ii).*gwy(jj);
            energynorm = energynorm + estress'*invC*estress*J.*scalefac.*gwx(ii).*gwy(jj);
            l2errvar = l2errvar + ((exdispmatx(jj,ii)-dispmatx(jj,ii))^2+(exdispmaty(jj,ii)-dispmaty(jj,ii))^2).*J.*scalefac.*gwx(ii).*gwy(jj);
            l2norm = l2norm + ((exdispmatx(jj,ii))^2+(exdispmaty(jj,ii))^2).*J.*scalefac.*gwx(ii).*gwy(jj);
            
        end
    end
    
    energyerrnorm = energyerrnorm + temp;

end

l2errvar = sqrt(l2errvar);
l2errrel = l2errvar/sqrt(l2norm);
energyerrnorm = sqrt(0.5*energyerrnorm);
EU = 0.5.*energynorm;
relenergynorm = energyerrnorm/sqrt(EU);


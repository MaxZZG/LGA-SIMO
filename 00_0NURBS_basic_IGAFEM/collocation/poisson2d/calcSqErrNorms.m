function [l2normerr,h1normerr,h2normerr]=calcSqErrNorms2(sol0, element_int, elementcounter, p, q, knotu, knotv, b_net, element_nod, coord_ij, nument)
% calculate the L2, H1 and H2 norms of the relative error for NURBS basis functions


num_gauss_x = p+1;
num_gauss_y = q+1;

[gwx, gpx]=quadrature(num_gauss_x, 'GAUSS', 1);
[gwy, gpy]=quadrature(num_gauss_y, 'GAUSS', 1);

l2norm = 0;
h1norm = 0;
h2norm = 0;
l2normerr = 0;
h1normerr = 0;
h2normerr = 0;

dim = 2;
%pre-compute the B-splines and 1st and 2nd derivatives

deriv_order = 2; %we need 2nd order derivatives
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
    index = 0;
    %for each gauss point
    for ii=1:num_gauss_x
        for  jj=1:num_gauss_y
            %find the knot span index in each direction
            ni = coord_ij(element_nod(i,end),1);
            nj = coord_ij(element_nod(i,end),2);
                        
            index = index + 1;                
            glob_index_u = pt_index_u(ni, ii);
            glob_index_v = pt_index_v(nj, jj);            
            
            %calculate the weights and control points corresponding to the
            %current element            
            wgts = reshape(b_net(ni-p:ni, nj-q:nj, dim+1), nument, 1);
            cpts = reshape(b_net(ni-p:ni, nj-q:nj, 1:2), nument, dim);
            
            M = M_arr_u(:,:,glob_index_u);
            P = M_arr_v(:,:,glob_index_v);                       
            
           [R, dR, ddR] = nurbshape2d5(M,P,p,q,wgts);
             
         
            %calculate the coordinates in the physical space
            coord = R*cpts;
         
            
            %calculate the Jacobian of the transformation
            dxdxi = dR*cpts;
            J = det(dxdxi);
            dxdxi = dxdxi';
            
            % Set up the second derivatives matrix and the matrix of squared first derivatives
            d2xdxi2 = ddR*cpts;
            
            dxdxi2 = [dxdxi(1,1)^2, dxdxi(1,1)*dxdxi(1,2), dxdxi(1,2)^2;...
                2*dxdxi(1,1)*dxdxi(2,1), dxdxi(1,1)*dxdxi(2,2)+dxdxi(1,2)*dxdxi(2,1), 2*dxdxi(1,2)*dxdxi(2,2);...
                dxdxi(2,1)^2, dxdxi(2,1)*dxdxi(2,2) dxdxi(2,2)^2];
            
            % Solve for first derivatives in global coordinates 
            dR = dxdxi'\dR;
                        
            % Solve for second derivatives in global coordinates 
            ddR = dxdxi2'\(ddR - d2xdxi2*dR);       
            
            xvals = coord(1);
            yvals = coord(2);
            %evaluate the exact solution and derivatives
            [u, dudx, dudy, ~, ~, ddudx, ddudxdy, ddudy ] = exact_sol(xvals,yvals);
            
            %evaluate the computed solution and derivatives
            cursol = sol0(element_nod(i,:));
            uh = R*cursol;
            
            duhdx = dR(1,:)*cursol;
            duhdy = dR(2,:)*cursol;
            dduhdx = ddR(1,:)*cursol;
            dduhdy = ddR(3,:)*cursol;
            dduhdxdy = ddR(2,:)*cursol;
            
            l2norm = l2norm + u^2*gwx(ii)*gwy(jj)*scalefac*J;
            h1norm = h1norm + (dudx^2+dudy^2)*gwx(ii)*gwy(jj)*scalefac*J;
            h2norm = h2norm + (ddudx^2 + ddudxdy^2 + ddudy^2)*gwx(ii)*gwy(jj)*scalefac*J;
            
            l2normerr = l2normerr + (u-uh)^2*gwx(ii)*gwy(jj)*scalefac*J;
            h1normerr = h1normerr + ((dudx-duhdx)^2+(dudy-duhdy)^2)*gwx(ii)*gwy(jj)*scalefac*J;
            h2normerr = h2normerr + ((ddudx-dduhdx)^2 + (ddudxdy-dduhdxdy)^2 + (ddudy-dduhdy)^2)*gwx(ii)*gwy(jj)*scalefac*J;
            
        end
    end
    
end

l2normerr = sqrt(l2normerr)/sqrt(l2norm);
h1normerr = sqrt(h1normerr)/sqrt(h1norm);
h2normerr = sqrt(h2normerr)/sqrt(h2norm);




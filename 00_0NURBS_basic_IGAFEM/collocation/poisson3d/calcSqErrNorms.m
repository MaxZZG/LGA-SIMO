function [l2normerr,h1normerr,h2normerr]=calcSqErrNorms(sol0, element_int, elementcounter, p, q, r, knotu, knotv, knotw, b_net, element_nod, coord_ijk)
% calculate the L2, H1 and H2 norms of the relative error for NURBS basis functions


num_gauss_x = p+1;
num_gauss_y = q+1;
num_gauss_z = r+1;

dim = 3; %we have 3 dimensions

[gwx, gpx]=quadrature(num_gauss_x, 'GAUSS', 1);
[gwy, gpy]=quadrature(num_gauss_y, 'GAUSS', 1);
[gwz, gpz]=quadrature(num_gauss_z, 'GAUSS', 1);

l2norm = 0;
h1norm = 0;
h2norm = 0;
l2normerr = 0;
h1normerr = 0;
h2normerr = 0;
 
%loop throught the knot spans in each direction and pre-compute
%Bsplines and connectivity arrays

deriv_order = 2; 
[ pt_index_u, M_arr_u ] = makeIndexBspline( knotu, p, deriv_order, gpx );
[ pt_index_v, M_arr_v ] = makeIndexBspline( knotv, q, deriv_order, gpy );
[ pt_index_w, M_arr_w ] = makeIndexBspline( knotw, r, deriv_order, gpz );

for i=1:elementcounter
   
    %for each element calculate the value of the basis
    %functions
   

    ximin = element_int(i,1);
    etamin = element_int(i,2);
    zetamin = element_int(i,3);
    ximax = element_int(i,4);
    etamax = element_int(i,5);         
    zetamax = element_int(i,6);

                              
    scalefac = (ximax - ximin)*(etamax - etamin)*(zetamax - zetamin)/8;
    
    %find the knot span index in each direction
    ni = coord_ijk(element_nod(i,end),1);
    nj = coord_ijk(element_nod(i,end),2);
    nk = coord_ijk(element_nod(i,end),3);

    %calculate the weights and control points corresponding to the
    %current element

    wgts= zeros((p+1)*(q+1)*(r+1), 1);
    cpts = zeros((p+1)*(q+1)*(r+1), dim);
    
    
    
    icount = 0;
    for kk = 1:r+1 
        for jj = 1:q+1
            for ii = 1:p+1
                icount = icount+1;
                wgts(icount) = b_net(ni+ii-p-1,nj+jj-q-1,nk+kk-r-1,dim+1);
                cpts(icount, 1) = b_net(ni+ii-p-1,nj+jj-q-1,nk+kk-r-1,1);
                cpts(icount, 2) = b_net(ni+ii-p-1,nj+jj-q-1,nk+kk-r-1,2);
                cpts(icount, 3) = b_net(ni+ii-p-1,nj+jj-q-1,nk+kk-r-1,3);
            end
        end
    end
    
    %for each gauss point
    for ii=1:num_gauss_x
        for jj=1:num_gauss_y
            for kk=1:num_gauss_z
                
                glob_index_u = pt_index_u(ni, ii);
                glob_index_v = pt_index_v(nj, jj);
                glob_index_w = pt_index_w(nk, kk);

                M = M_arr_u(:,:,glob_index_u);
                P = M_arr_v(:,:,glob_index_v);
                Q = M_arr_w(:,:,glob_index_w);
                
                %compute the value of the shape function and gradient with respect to parameter space                         
                 [R, dR, ddR] = nurbshape3d5(M,P,Q,p,q,r,wgts);
                                
                % calculate coordinates in physical space
                coord = R*cpts;

                % Compute Jacobian matrix
                dxdxi = dR*cpts;
                dxdxi = dxdxi';
                J = det(dxdxi);
                

                % Set up the second derivatives matrix and the matrix of squared first derivatives
                d2xdxi2 = ddR*cpts;

                dxdxi2 = [ dxdxi(1,1)^2            dxdxi(1,2)^2             dxdxi(1,3)^2            dxdxi(1,1)*dxdxi(1,2)                       dxdxi(1,1)*dxdxi(1,3)                       dxdxi(1,2)*dxdxi(1,3)
                           dxdxi(2,1)^2            dxdxi(2,2)^2             dxdxi(2,3)^2            dxdxi(2,1)*dxdxi(2,2)                       dxdxi(2,1)*dxdxi(2,3)                       dxdxi(2,2)*dxdxi(2,3)
                           dxdxi(3,1)^2            dxdxi(3,2)^2             dxdxi(3,3)^2            dxdxi(3,1)*dxdxi(3,2)                       dxdxi(3,1)*dxdxi(3,3)                       dxdxi(3,2)*dxdxi(3,3)
                           2*dxdxi(1,1)*dxdxi(2,1) 2*dxdxi(1,2)*dxdxi(2,2)  2*dxdxi(1,3)*dxdxi(2,3) dxdxi(1,1)*dxdxi(2,2)+dxdxi(1,2)*dxdxi(2,1) dxdxi(1,1)*dxdxi(2,3)+dxdxi(1,3)*dxdxi(2,1) dxdxi(1,2)*dxdxi(2,3)+dxdxi(1,3)*dxdxi(2,2)
                           2*dxdxi(1,1)*dxdxi(3,1) 2*dxdxi(1,2)*dxdxi(3,2)  2*dxdxi(1,3)*dxdxi(3,3) dxdxi(1,1)*dxdxi(3,2)+dxdxi(1,2)*dxdxi(3,1) dxdxi(1,1)*dxdxi(3,3)+dxdxi(1,3)*dxdxi(3,1) dxdxi(1,2)*dxdxi(3,3)+dxdxi(1,3)*dxdxi(3,2)
                           2*dxdxi(2,1)*dxdxi(3,1) 2*dxdxi(2,2)*dxdxi(3,2)  2*dxdxi(2,3)*dxdxi(3,3) dxdxi(2,1)*dxdxi(3,2)+dxdxi(2,2)*dxdxi(3,1) dxdxi(2,1)*dxdxi(3,3)+dxdxi(2,3)*dxdxi(3,1) dxdxi(2,2)*dxdxi(3,3)+dxdxi(2,3)*dxdxi(3,2) ];

                % Solve for first derivatives in global coordinates 
                dR = dxdxi'\dR;

                % Solve for second derivatives in global coordinates 
                ddR = dxdxi2'\(ddR - d2xdxi2*dR);

                xvals = coord(1);
                yvals = coord(2);             
                zvals = coord(3);
                %evaluate the exact solution and derivatives
                [ u, dudx, dudy, dudz, ~, ~, ddudx, ddudy, ddudz, ddudxdy, ddudxdz, ddudydz ] = exact_sol(xvals,yvals,zvals);

                %evaluate the computed solution and derivatives
                cursol = sol0(element_nod(i,:));                                                            
                uh = R*cursol;

                duhdx = dR(1,:)*cursol;
                duhdy = dR(2,:)*cursol;
                duhdz = dR(3,:)*cursol;
                dduhdx = ddR(1,:)*cursol;
                dduhdy = ddR(2,:)*cursol;
                dduhdz = ddR(3,:)*cursol;
                dduhdxdy = ddR(4,:)*cursol;                                            
                dduhdxdz = ddR(5,:)*cursol;                                            
                dduhdydz = ddR(6,:)*cursol;                                            

                l2norm = l2norm + u^2*gwx(ii)*gwy(jj)*gwz(kk)*scalefac*J;
                h1norm = h1norm + (dudx^2+dudy^2+dudz^2)*gwx(ii)*gwy(jj)*gwz(kk)*scalefac*J;
                h2norm = h2norm + (ddudx^2+ ddudy^2 + ddudz^2 + ddudxdy^2 + ddudxdz^2 + ddudydz^2 )...
                    *gwx(ii)*gwy(jj)*gwz(kk)*scalefac*J;

                l2normerr = l2normerr + (u-uh)^2*gwx(ii)*gwy(jj)*gwz(kk)*scalefac*J;
                h1normerr = h1normerr + ((dudx-duhdx)^2+(dudy-duhdy)^2+(dudz-duhdz)^2)*gwx(ii)*gwy(jj)*gwz(kk)*scalefac*J;
                h2normerr = h2normerr + ((ddudx-dduhdx)^2 + (ddudy-dduhdy)^2 + (ddudz-dduhdz)^2 + (ddudxdy-dduhdxdy)^2 +...
                    (ddudxdz-dduhdxdz)^2 + (ddudydz-dduhdydz)^2)*gwx(ii)*gwy(jj)*gwz(kk)*scalefac*J;
            end
        end
    end                                                   
    
end

l2normerr = sqrt(l2normerr)/sqrt(l2norm);
h1normerr = sqrt(h1normerr)/sqrt(h1norm);
h2normerr = sqrt(h2normerr)/sqrt(h2norm);

    
    
    
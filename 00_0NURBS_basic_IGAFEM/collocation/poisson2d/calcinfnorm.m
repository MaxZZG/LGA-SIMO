function [errw0inf, errw1inf, errw2inf] = calcinfnorm2(sol0, elementcounter, p, q, knotu, knotv, b_net, element_nod, coord_ij, nument)
%calculates the infinity norm of the error in the solution

numpts = 10; %number of sample points per element

errw0inf = 0;
errw1inf = 0;
errw2inf = 0;
w0inf = 0;
w1inf = 0;
w2inf = 0;

dim = 2;

%pre-compute the B-splines and 1st and 2nd derivatives

px = linspace(-1, 1, numpts);
py = linspace(-1, 1, numpts);

deriv_order = 2; %we need 2nd order derivatives
[ pt_index_u, M_arr_u ] = makeIndexBspline( knotu, p, deriv_order, px );
[ pt_index_v, M_arr_v ] = makeIndexBspline( knotv, q, deriv_order, py );


for i=1:elementcounter
    
    
    uh = zeros(numpts, numpts);
    duhdx = zeros(numpts, numpts);
    duhdy = zeros(numpts, numpts);
    dduhdx = zeros(numpts, numpts);
    dduhdy = zeros(numpts, numpts);
    dduhdxdy = zeros(numpts, numpts);
    %  uex = zeros(numpts, numpts);
    %  duhdxex = zeros(numpts, numpts);
    xvals = zeros(numpts, numpts);
    yvals = zeros(numpts, numpts);
    
    index = 0;
    %for each plot point
    for ii=1:numpts
        for jj=1:numpts                                    
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
            
            %Strore the values of the coordinates
            xvals(jj,ii) = coord(1);
            yvals(jj,ii) = coord(2);
            
            cursol = sol0(element_nod(i,:));
            uh(jj,ii) = R*cursol;
            duhdx(jj,ii) = dR(1,:)*cursol;
            duhdy(jj,ii) = dR(2,:)*cursol;
            dduhdx(jj,ii) = ddR(1,:)*cursol;
            dduhdy(jj,ii) = ddR(3,:)*cursol;
            dduhdxdy(jj,ii) = ddR(2,:)*cursol;
        end
    end
    
    
    %evaluate the exact solution
    [u, dudx, dudy, ~, ~, ddudx, ddudxdy, ddudy ] = exact_sol(xvals,yvals);
    
    w0inf = max([w0inf, max(max(abs(u)))]);
    w1inf = max([w1inf, max(max(abs(dudx))), max(max(abs(dudy)))]);
    w2inf = max([w2inf, max(max(abs(ddudx))), max(max(abs(ddudxdy))), max(max(abs(ddudy)))]);
    
    errw0inf = max([errw0inf, max(max(abs(u-uh)))]);
    errw1inf = max([errw1inf, max(max(abs(dudx - duhdx))), max(max(abs(dudy - duhdy)))]);
    errw2inf = max([errw2inf, max(max(abs(ddudx-dduhdx))), max(max(abs(ddudxdy-dduhdxdy))), max(max(abs(ddudy-dduhdy)))]);
end

errw0inf = errw0inf/w0inf;
errw1inf = errw1inf/w1inf;
errw2inf = errw2inf/w2inf;



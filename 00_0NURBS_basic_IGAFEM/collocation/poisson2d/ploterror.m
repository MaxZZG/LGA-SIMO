function ploterror(sol0, elementcounter, p, q, knotu, knotv, b_net, lenu, lenv, element_nod, coord_ij)
%plots the error in the solution

numpts = 5; %number of plot points per element
nument = (p+1)*(q+1);

for i=1:elementcounter
    
    
    %for each integration element calculate the value of the basis
    %functions
    
    
    px = linspace(-1, 1, numpts);
    py = linspace(-1, 1, numpts);
    
    uh = zeros(numpts, numpts);
    duhdx = zeros(numpts, numpts);
    duhdy = zeros(numpts, numpts);
    dduhddx = zeros(numpts, numpts);
    dduhddy = zeros(numpts, numpts);
    dduhdxdy = zeros(numpts, numpts);
    xvals = zeros(numpts, numpts);
    yvals = zeros(numpts, numpts);
    
    
    %for each plot point
    for ii=1:numpts
        for jj=1:numpts
            %compute the value of the shape function and gradient
            [R, dR, ddR, coord] = nurbshaped(i, px(ii), py(jj), knotu, knotv, b_net, p, q, lenu, lenv, element_nod, coord_ij);
            xvals(jj,ii) = coord(1);
            yvals(jj,ii) = coord(2);
            
            for j=1:nument
                cursol = sol0(element_nod(i,j));
                uh(jj,ii) = uh(jj,ii) + cursol*R(j);
                duhdx(jj,ii) = duhdx(jj,ii) + cursol*dR(1,j);
                duhdy(jj,ii) = duhdy(jj,ii) + cursol*dR(2,j);
                dduhddx(jj,ii) = dduhddx(jj,ii) + cursol*ddR(1,j);
                dduhddy(jj,ii) = dduhddy(jj,ii) + cursol*ddR(3,j);
                dduhdxdy(jj,ii) = dduhdxdy(jj,ii) + cursol*ddR(2,j);
                
            end
        end
    end
    
    %evaluate the exact solution
    [u, dudx, dudy, f, alpha, ddudx, ddudxdy, ddudy] = exact_sol(xvals,yvals);
    surf(xvals,yvals,u-uh)
    
    %pause
    hold on
    
end
title('u-u_h')
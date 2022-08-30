function ploterr(sol0, span_elm, p, q, r, knotu, knotv, knotw, b_net, element_nod, coord_ijk, lenu, lenv, lenw)

figure

dim = 3; %we have 3 dimensions

w_val = 0.5;
numpts = 100;
fudge = 1e-4;
px = linspace(fudge,1-fudge,numpts);
py = linspace(fudge,1-fudge,numpts);
%pz = linspace(fudge,1-fudge,numpts);

num_plot_pts = numpts^2;


%establish sample point -> element connectivity
plot_elm = zeros(num_plot_pts, 1);
plot_pts = zeros(num_plot_pts, dim);
plot_index = zeros(num_plot_pts, 2); %store the ii,jj index for plotting
index = 0;
for i=1:numpts
    for j=1:numpts
        index = index + 1;
        curu = px(i);
        curv = py(j);
        curw = w_val;
        plot_pts(index,:) = [curu, curv, curw];
        plot_index(index, :) = [i, j];
        for t1 = 1:length(knotu)-1
            if (curu >= knotu(t1)) && (curu <= knotu(t1+1)) && (knotu(t1) < knotu(t1+1))
                for t2 = 1:length(knotv)-1
                    if (curv >= knotv(t2)) && (curv <= knotv(t2+1)) && (knotv(t2) < knotv(t2+1))     
                        for t3 = 1:length(knotw)-1
                            if (curw >= knotw(t3)) && (curw <= knotw(t3+1)) && (knotw(t3) < knotw(t3+1))
                                plot_elm(index) = span_elm(t1,t2,t3);                  
                            end
                        end
                    end
                end
            end
        end
    end
end

uh = zeros(numpts, numpts);
duhdx = zeros(numpts, numpts);
duhdy = zeros(numpts, numpts);     
duhdz = zeros(numpts, numpts);

uex = zeros(numpts,numpts);
dudx = zeros(numpts, numpts);
dudy = zeros(numpts, numpts);
dudz = zeros(numpts, numpts);

xvals = zeros(numpts, numpts);
yvals = zeros(numpts, numpts);
zvals = zeros(numpts, numpts);

for i_plt=1:num_plot_pts
   
    %for plot point calculate the value of the basis
    %functions
    i = plot_elm(i_plt);      
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
    
    ii = plot_index(i_plt,1);
    jj = plot_index(i_plt,2);


    u = plot_pts(i_plt,1);
    v = plot_pts(i_plt,2);
    w = plot_pts(i_plt,3);
    

    %compute the value of the shape function and gradient with respect to parameter space                         
    [R, dRdx, ~, coord, J] = nurbshape3da(i, u, v, w, knotu, knotv, knotw, b_net, p, q, r, lenu, lenv, lenw, element_nod, coord_ijk);

    xvals(jj,ii) = coord(1);
    yvals(jj,ii) = coord(2);             
    zvals(jj,ii) = coord(3);
    %evaluate the exact solution and derivatives
    [ uex(jj,ii), dudx(jj,ii), dudy(jj,ii), dudz(jj,ii)] = exact_sol(coord(1),coord(2),coord(3));


    %evaluate the computed solution and derivatives
    cursol = sol0(element_nod(i,:));                                                            
    uh(jj,ii) = R*cursol;

   
   % u-uh
   % pause

    %pause
    hold on
                         
    
end


surf(xvals,yvals,zvals, (uex-uh)/max(max(abs(uex))))
colorbar
shading interp
view([80, 20])
title('Error (u-u_h)/|u|')

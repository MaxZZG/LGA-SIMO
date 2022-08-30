%script for plotting the solution and the errors for the "Plate with a
%hole" problem

%Displaying the displacements
element4 = zeros(elementcounter, 4);
physcoord = zeros(4*elementcounter, 2);
dispcoord = zeros(4*elementcounter, 2);
straincoord = zeros(4*elementcounter, 3);
sigmacoord = zeros(4*elementcounter, 7);


numix = length(knotu) - 1 - 2*p;
numiy = length(knotv) - 1 - 2*q;

streest = zeros(numix, numiy);

for i=1:elementcounter
    
    element4(i, :) = [(i-1)*4+1:(i-1)*4+4];    
    
    nx = 2;
    ny = 2;
    px = linspace(-1, 1, nx);
    py = linspace(-1, 1, ny);
   
    
    ximin = element_int(i,1);
    etamin = element_int(i,2);
    ximax = element_int(i,3);
    etamax = element_int(i,4);
    
    coordt = cell(ny,nx);
    
    dispmatx = zeros(ny,nx);
    dispmaty = zeros(ny,nx);
    strain11 = zeros(ny,nx);
    strain12 = zeros(ny,nx);
    strain22 = zeros(ny,nx);
    stressvect = cell(ny,nx);
    exstressvect = cell(ny,nx);
    
    sigmad = zeros(ny,nx);
    
    scrtx = element_nod(i,:);
    %reverse order of scrtx to match format of nurb
    scrtx = scrtx(end:-1:1);
    
    for ii=1:nx
        for jj=1:ny
            
            upt =((ximax-ximin)*px(ii) + ximax + ximin)/2;
            vpt =((etamax-etamin)*py(jj) + etamax + etamin)/2;
            [R, coord, dRdxi, dRdx, J] = nurbshape(i, upt, vpt, knotu, knotv, b_net, p, q, lenu, lenv, element_nod, coord_ij);
            
            coordt{jj,ii} = coord;
            
            cdRdx = dRdx(:,1);
            cdRdy = dRdx(:,2);
                        
            %calculate displacement values
            dispmatx(jj,ii) = dispmatx(jj,ii) + R'*x(2*scrtx-1);
            dispmaty(jj,ii) = dispmaty(jj,ii) + R'*x(2*scrtx);

            %calculate strain values
            strain11(jj,ii) = strain11(jj,ii) + cdRdx'*x(2*scrtx-1);
            strain12(jj,ii) = strain12(jj,ii) + cdRdy'*x(2*scrtx-1) + cdRdx'*x(2*scrtx);
            strain22(jj,ii) = strain22(jj,ii) + cdRdy'*x(2*scrtx);
      
            stressvect{jj,ii} = C*[strain11(jj,ii); strain22(jj,ii); strain12(jj,ii)];
            
            %calculate the exact stresses
            physx = coord(1);
            physy = coord(2);
            exstressvect{jj,ii} = ghole(physx, physy, rad, tx);
                        
            sigma11 = stressvect{jj,ii}(1);
            sigma12 = stressvect{jj,ii}(3);
            sigma22 = stressvect{jj,ii}(2);
            
            lambda = Emod*nu/(1-nu^2);
            mu = Emod/(2*(1+nu));
            
            sigmad(jj,ii) = 1/(4*mu)*((mu^2/(6*(mu+lambda)^2)+1/2)*(sigma11+sigma22)^2+2*(sigma12^2-sigma11*sigma22));
            
        end
    end
        
    
    physcoord((i-1)*4+1, :) = coordt{1,1};
    physcoord((i-1)*4+2, :) = coordt{1,2};
    physcoord((i-1)*4+3, :) = coordt{2,2};
    physcoord((i-1)*4+4, :) = coordt{2,1};
    
    icoord = coord_ij(element_nod(i,1),1);
    jcoord = coord_ij(element_nod(i,1),2);
    
    
    maxsxx = max([stressvect{1,1}(1), stressvect{1,2}(1), stressvect{2,1}(1), stressvect{2,2}(1)]);
    minsxx = min([stressvect{1,1}(1), stressvect{1,2}(1), stressvect{2,1}(1), stressvect{2,2}(1)]);
    
    maxsyy = max([stressvect{1,1}(2), stressvect{1,2}(2), stressvect{2,1}(2), stressvect{2,2}(2)]);
    minsyy = min([stressvect{1,1}(2), stressvect{1,2}(2), stressvect{2,1}(2), stressvect{2,2}(2)]);
    
    maxsxy = max([stressvect{1,1}(3), stressvect{1,2}(3), stressvect{2,1}(3), stressvect{2,2}(3)]);
    minsxy = min([stressvect{1,1}(3), stressvect{1,2}(3), stressvect{2,1}(3), stressvect{2,2}(3)]);
    
    
    streest(icoord, jcoord) = (maxsxx - minsxx)^2 + (maxsyy - minsyy)^2 + (maxsxy - minsxy)^2;
    
    
    dispcoord((i-1)*4+1, :) = [dispmatx(1,1) dispmaty(1,1)];
    dispcoord((i-1)*4+2, :) = [dispmatx(1,2) dispmaty(1,2)];
    dispcoord((i-1)*4+3, :) = [dispmatx(2,2) dispmaty(2,2)];
    dispcoord((i-1)*4+4, :) = [dispmatx(2,1) dispmaty(2,1)];
    straincoord((i-1)*4+1, :) = [strain11(1,1) strain12(1,1) strain22(1,1)];
    straincoord((i-1)*4+2, :) = [strain11(1,2) strain12(1,2) strain22(1,2)];
    straincoord((i-1)*4+3, :) = [strain11(2,2) strain12(2,2) strain22(2,2)];
    straincoord((i-1)*4+4, :) = [strain11(2,1) strain12(2,1) strain22(2,1)];
    
    sigmacoord((i-1)*4+1, :) = [stressvect{1,1}' sigmad(1,1) exstressvect{1,1}'];
    sigmacoord((i-1)*4+2, :) = [stressvect{1,2}' sigmad(1,2) exstressvect{1,2}'];
    sigmacoord((i-1)*4+3, :) = [stressvect{2,2}' sigmad(2,2) exstressvect{2,2}'];
    sigmacoord((i-1)*4+4, :) = [stressvect{2,1}' sigmad(2,1) exstressvect{2,1}'];
end



%maginifcation factor for the displacement plot
factor = 1000;

%plot the computed stresses

figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,1), 'facecolor','interp', 'EdgeColor','none')
view(0,90)
title('Displacements and \sigma_{11}')
colorbar('vert')

figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,3), 'facecolor','interp', 'EdgeColor','none')
view(0,90)
title('Displacements and \sigma_{12}')
colorbar('vert')

figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,2), 'facecolor','interp', 'EdgeColor','none')
view(0,90)
title('Displacements and \sigma_{22}')
colorbar('vert')

%plot the exact stresses

figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,5), 'facecolor','interp', 'EdgeColor','none')
view(0,90)
title('Displacements and exact \sigma_{11}')
colorbar('vert')

figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,7), 'facecolor','interp', 'EdgeColor','none')
view(0,90)
title('Displacements and exact \sigma_{12}')
colorbar('vert')

figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,6), 'facecolor','interp', 'EdgeColor','none')
view(0,90)
title('Displacements and exact \sigma_{22}')
colorbar('vert')


%plot the errors
figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,5)-sigmacoord(:,1), 'facecolor','interp', 'EdgeColor','none')
view(0,90)
title('Displacements and error in \sigma_{11}')
colorbar('vert')

figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,7)-sigmacoord(:,3), 'facecolor','interp', 'EdgeColor','none')
view(0,90)
title('Displacements and error in \sigma_{12}')
colorbar('vert')

figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,6)-sigmacoord(:,2), 'facecolor','interp', 'EdgeColor','none')
view(0,90)
title('Displacements and error in \sigma_{22}')
colorbar('vert')


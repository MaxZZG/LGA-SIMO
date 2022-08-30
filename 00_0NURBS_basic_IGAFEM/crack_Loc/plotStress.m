% compute stresses and displacements at nodes of
% the visualization mesh
% Also export this mesh together with stresses
% and displacements to VTK file which can then be processed by Paraview.
% Vinh Phu Nguyen
% Delft University of Technology

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute stresses and displacements

stress = zeros(noElems,size(elementV,2),4);
displacement   = zeros(noElems,size(elementV,2),2);

for e=1:numelem
    ni = inn(ien(e,1),1);% get NURBS coordinates
    nj = inn(ien(e,1),2);
    if (abs(uKnot(ni)-uKnot(ni+1))>tol)&&(abs(vKnot(nj)-vKnot(nj+1))>tol)
        xiE=[uKnot(ni), uKnot(ni+1)];
        etaE=[vKnot(nj), vKnot(nj+1)];
        
        sctr   = element(e,:);          %  element scatter vector
        nn     = length(sctr);
        pts    = controlPts(sctr,:);
    
        uspan = FindSpan(noPtsX-1,p,xiE(1),uKnot);
        vspan = FindSpan(noPtsY-1,q,etaE(1),vKnot);
    
        elemDisp  = element_disp(e,pos,enrich_node,U);
    
        
%     idu    = index(e,1);
%     idv    = index(e,2);
%     xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
%     etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
%     
%     sctr   = element(e,:);          %  element scatter vector
%     nn     = length(sctr);
%     
%     pts    = controlPts(sctr,:);
%     
%     uspan = FindSpan(noPtsX-1,p,xiE(1),uKnot);
%     vspan = FindSpan(noPtsY-1,q,etaE(1),vKnot);
%     
%     elemDisp  = element_disp(e,pos,enrich_node,U);
    
    % loop over Gauss points
    
    gp = 1;
    for iv=1:2
        if (iv==2)
            xiE = sort(xiE,'descend');
        end
        for iu=1:2            
            Xi  = xiE(iu);
            Eta = etaE(iv);
            
            [N dRdxi dRdeta] = NURBS2DBasisDersSpecial([Xi; Eta],...
                p,q,uKnot,vKnot,weights',[uspan;vspan]);
            
            [B,w] = BMatrixXIGA(Xi,Eta,e,enrich_node,N, dRdxi, dRdeta);
            [exN] = NMatrixXIGA(e,enrich_node,N);
            
            strain          = B*elemDisp;
            sigma           = C*strain;
            stress(e,gp,1:3)= sigma;
            
            % von Mises stress
            stress(e,gp,4)  = sqrt(sigma(1)^2+sigma(2)^2-...
                                   sigma(1)*sigma(2)+3*sigma(3)^2); 
            
            % the following is incorrect since
            % enriched dofs are not included!!!
            %disp(e,gp,:)    = N*[Ux(sctr) Uy(sctr)];
            
            displacement(e,gp,:)    = exN*[elemDisp(1:2:end) ...
                                   elemDisp(2:2:end)];            
            gp = gp +1;
        end
    end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% export to VTK format to plot in Mayavi or Paraview

sigmaXX = zeros(size(node,1),1);
sigmaYY = zeros(size(node,1),1);
sigmaXY = zeros(size(node,1),1);
sigmaVM = zeros(size(node,1),1);

dispX = zeros(size(node,1),1);
dispY = zeros(size(node,1),1);

for e=1:size(elementV,1)
    connect = elementV(e,:);
    for in=1:4
        nid = connect(in);
        sigmaXX(nid) = stress(e,in,1);
        sigmaYY(nid) = stress(e,in,2);
        sigmaXY(nid) = stress(e,in,3);
        sigmaVM(nid) = stress(e,in,4);
        
        dispX(nid) = displacement(e,in,1);
        dispY(nid) = displacement(e,in,2);
    end
end

VTKPostProcess(node,elementV,'Quad4',vtuFile,...
               [sigmaXX sigmaYY sigmaXY sigmaVM],[dispX dispY]);


stressComp=2;
figure
clf
hold on
% plot_field(node+30*[dispX dispY],elementV,'Q4',sigmaYY);  % ve voi bien dang
plot_field(node+node,elementV,'Q4',sigmaYY);  % ve voi hinh dang ban dau
hold on
colorbar
title('Stress in x direction')
axis off
%plot_mesh(node,elementV,'Q4','k.-');

figure
clf
hold on
% plot_field(node+30*[dispX dispY],elementV,'Q4',sigmaVM);
plot_field(node+node,elementV,'Q4',sigmaVM);
hold on
colorbar
title('Von-Mises Stress')
axis off
% plot_mesh(node,elementV,'Q4','k.-');

figure
clf
plot_field(node,elementV,'Q4',displacement(:,:,2));
hold on
colorbar
title('Displacement in y direction')
axis off

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)


figure
clf
plot_field(node,elementV,'Q4',displacement(:,:,1));
hold on
colorbar
title('Displacement in x direction')
axis off

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)

clear disp;

%
           
 


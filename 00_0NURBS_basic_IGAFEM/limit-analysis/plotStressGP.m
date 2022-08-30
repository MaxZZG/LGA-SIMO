% compute stresses and displacements at nodes of
% the visualization mesh
% Also export this mesh together with stresses
% and displacements to VTK file which can then be processed by Paraview.
% Vinh Phu Nguyen
% Delft University of Technology

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute stresses and displacements

% stress = zeros(noElems,size(elementV,2),4);
% disp   = zeros(noElems,size(elementV,2),2);

count=1; coordGP=[];
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
        
        
    if     (ismember(e,split_elem)) % split element, 13 GPs/subtriangle
%         [W,Q] = discontQ4quad(7,levelSets(1,elementV(e,:),1));
        [W,Q] = quadrature(20,'GAUSS',2);
    elseif (ismember(e,tip_elem))   % tip element
%         [W,Q] = disTipQ4quad(7,levelSets(1,elementV(e,:)),node(sctrV,:),xTip);
        [W,Q] = quadrature(20,'GAUSS',2);
    elseif ( any(intersect(tip_nodes,sctr)) ~= 0)% having tip enriched nodes
        [W,Q] = quadrature(10,'GAUSS',2);
    else
        [W,Q] = quadrature(noGPs1,'GAUSS',2);
    end
    
    
    % loop over Gauss points
         for gp=1:size(W,1)
            pt      = Q(gp,:);
        % compute coords in parameter space
            Xi   = ((xiE(2) - xiE(1))*pt(1)+xiE(2) + xiE(1))/2;
            Eta  = ((etaE(2)-etaE(1))*pt(2)+etaE(2)+etaE(1))/2;
            [N dRdxi dRdeta] = NURBS2DBasisDersSpecial([Xi; Eta],...
                               p,q,uKnot,vKnot,weights',[uspan;vspan]);
            
            [B,w] = BMatrixXIGA(Xi,Eta,e,enrich_node,N, dRdxi, dRdeta);
            [exN] = NMatrixXIGA(e,enrich_node,N);
            
            strain          = B*elemDisp;
            sigma           = C*strain;
            stress(count,1:3)= sigma;
            
            % von Mises stress
            stress(count,4)  = sqrt(sigma(1)^2+sigma(2)^2-...
                                   sigma(1)*sigma(2)+3*sigma(3)^2); 
            
            % the following is incorrect since
            % enriched dofs are not included!!!
            %disp(e,gp,:)    = N*[Ux(sctr) Uy(sctr)];
            
            disp(count,:)    = exN*[elemDisp(1:2:end) ...
                                   elemDisp(2:2:end)];            
            Gpt    = N*controlPts(sctr,:);
            
            coordGP= [coordGP;Gpt];
            count=count+1;
         end
    end
end

% get decompused triangles
tri = delaunay(coordGP(:,1),coordGP(:,2));
tri = tricheck(coordGP,tri);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% export to VTK format to plot in Mayavi or Paraview

sigmaXX = stress(:,1);
sigmaYY = stress(:,2);
sigmaXY = stress(:,3);
sigmaVM = stress(:,4);

dispX = disp(:,1);
dispY = disp(:,2);


clear stress disp

VTKPostProcess(coordGP,tri,2,'Tri3',vtuFile,...
               [sigmaXX sigmaYY sigmaXY sigmaVM],[dispX dispY]);           

% plot stress

stressComp=2;
figure
clf
hold on
plot_field(coordGP+30*[dispX dispY],tri,'T3',sigmaYY);   % stress+deformation
hold on
colorbar
title('Stress in x direction')
axis off; view(2);
% plot_mesh(node,elementV,'Q4','k.-');

figure
clf
hold on
plot_field(coordGP+30*[dispX dispY],tri,'T3',sigmaVM);   % stress+deformation
hold on
colorbar
title('Von-Mises Stress')
axis off; view(2);
% plot_mesh(node,elementV,'Q4','k.-');


figure
clf
plot_field(coordGP,tri,'T3',dispY);
hold on
colorbar
title('Displacement in y direction')
axis off; view(2);

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)


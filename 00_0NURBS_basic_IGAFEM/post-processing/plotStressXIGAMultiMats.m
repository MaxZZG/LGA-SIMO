% compute 2D stresses and displacements at nodes of
% the visualization mesh
% Also export this mesh together with stresses/strains
% and displacements to VTK file which can then be processed by Paraview.
% This is for problems with multiple materials
% Vinh Phu Nguyen
% Delft University of Technology

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute stresses and displacements

stress = zeros(noElems,4,4);
strain = zeros(noElems,4,3);
displacement   = zeros(noElems,4,2);

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);          %  element scatter vector of NURBS mesh
    sctrV  = elementV(e,:);         %  element scatter vector of Q4 FE mesh
    nn     = length(sctr);
    
    pts    = controlPts(sctr,:);
    
    uspan = FindSpan(noPtsX-1,p,xiE(1), uKnot);
    vspan = FindSpan(noPtsY-1,q,etaE(1),vKnot);
    
    elemDisp  = element_disp(e,pos,enrich_node,U);
    
    % loop over Gauss points
    
    gp = 1;
    for iv=1:2
        if (iv==2)
            xiE = sort(xiE,'descend'); % counter-clockwise node numbering
        end
        for iu=1:2
            Xi  = xiE(iu);
            Eta = etaE(iv);
            
            [N, dRdxi, dRdeta] = NURBS2DBasisDersSpecial([Xi; Eta],...
                p,q,uKnot,vKnot,weights',[uspan;vspan]);
            
            [B,w] = BMatrixXIGA(Xi,Eta,e,enrich_node,N,dRdxi,dRdeta);
            [exN] = NMatrixXIGA(e,enrich_node,N);
            
            %if (node(sctrV(gp),2) - L/2 >= 0)
            if (chi(sctrV(gp)) >= 0)
                C = Cm;
            else
                C = Ci;
            end
            
            strainGp        = B*elemDisp;
            sigma           = C*strainGp;
            stress(e,gp,1:3)= sigma;
            strain(e,gp,1:3)= strainGp;
            
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% export to VTK format to plot in Mayavi or Paraview
noNodes = size(node,1);

sigmaXX = zeros(noNodes,1);
sigmaYY = zeros(noNodes,1);
sigmaXY = zeros(noNodes,1);
sigmaVM = zeros(noNodes,1);

epsXX = zeros(noNodes,1);
epsYY = zeros(noNodes,1);
epsXY = zeros(noNodes,1);

dispX = zeros(noNodes,1);
dispY = zeros(noNodes,1);

for e=1:size(elementV,1)
    connect = elementV(e,:);
    for in=1:4
        nid          = connect(in);
        sigmaXX(nid) = stress(e,in,1);
        sigmaYY(nid) = stress(e,in,2);
        sigmaXY(nid) = stress(e,in,3);
        sigmaVM(nid) = stress(e,in,4);
        
        epsXX(nid) = strain(e,in,1);
        epsYY(nid) = strain(e,in,2);
        epsXY(nid) = strain(e,in,3);
        
        dispX(nid) = displacement(e,in,1);
        dispY(nid) = displacement(e,in,2);
    end
end

VTKPostProcess(node,elementV,'Quad4',vtuFile,...
    [sigmaXX sigmaYY sigmaXY sigmaVM],[dispX dispY]);




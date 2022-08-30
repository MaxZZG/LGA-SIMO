% data for edge crack problem with C2 elements
p = 2;
q = p;
cont=p-1;   % continuity derivative
Rep=p-cont; % repeated knot
% 
Numy = 12;   % numr 
Numx = 60;   % numtheta,  luon la so le de vet nut cat qua phtu!!!!!
% 
R_out = 60 ;
r_in  = 50 ;
theta=90*pi/180;
%%%%%%%% section annular plate
[CP,U,V,p,q]=pipe_coasemesh(R_out,r_in,p,theta);
R1 = refinement_vec_repeated(U,Numx,p-cont);
R2 = refinement_vec_repeated(V,Numy,q-cont);

% R1 = refinement_vec(U,ref);
% R2 = refinement_vec(V,ref);

[CP,uKnot,vKnot] = knot_refine_surf(p,q,U,V,CP,R1,R2);

plotNURBS_surf_El_CP(p,q,uKnot,vKnot,CP); hold on
view(2)

% plotNURBS_surf_El_CP(p,q,uKnot,vKnot,CP); hold on
plot_ctrlnet(CP,'ro');view(2)

noPtsX = size(CP,1);
noPtsY = size(CP,2);
noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

weights = reshape(CP(:,:,4),noCtrPts,1);

gcoord(:,1)=reshape(CP(:,:,1),noCtrPts,1);
gcoord(:,2)=reshape(CP(:,:,2),noCtrPts,1);
controlPts=gcoord;

% generate connectivities and index for control point   
[ien,inn]=genIEN_INN_2D_repeated_Ck(p,q,noPtsX,noPtsY);
element=sort(ien,2);


% generate element connectivity ...
[ien,inn]=genIEN_INN_2D_repeated_Ck(p,q,noPtsX,noPtsY);
element=sort(ien,2);

% generate connectivities and index for physis element
buildVisualizationMesh;

% crack data

noCracks = 1;                 % number of cracks

xCrack   = zeros(noCracks,2,2);
xTips    = zeros(noCracks,2);

a = (R_out-r_in)/2;                     % crack length

xCr   = [r_in 0; r_in+a 0];
xTip  = [r_in+a 0];
seg   = xCr(2,:) - xCr(1,:);   % tip segment
alpha = atan2(seg(2),seg(1));  % inclination angle
QT    =[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];

xCrack(1,:,:) = xCr;
xTips(1,:)    = xTip; 


% level set computation

numnode   = size(node,1);
numelem   = size(elementV,1);

levelSetCracks

% Choose enriched nodes...
% chooseEnrichedNodes

%
enrich_node = zeros(noCtrPts,1);
crack_node  = zeros(noCtrPts,1); % which crack to which the node enriched

count1 = 0;
count2 = 0;
tol = 1e-8;

tNodes = 1;
for i=1:p
    tNodes = [tNodes (p+1)*i+1]; 
end

% loop over elements
for iel = 1 : numelem
    sctr    = elementV(iel,:);
    sctrIGA = element(iel,:);
    
    % loop over cracks
    for iCr = 1 : noCracks
        phi  = levelSets(iCr,sctr,1);
        psi  = levelSets(iCr,sctr,2);
        
        if ( max(phi)*min(phi) <= tol ) && (max(psi)*min(psi) < tol)
            count2               = count2 + 1 ; % ah, one tip element
            tip_elem(count2)     = iel;
            enrich_node(sctrIGA) = 2;
            crack_node(sctrIGA)  = iCr;
        end
                
        if ( max(phi)*min(phi) < tol ) 
            if (max(psi) < 0) && (min(node(sctr,1)) > 0)
                count1                 = count1 + 1 ; % ah, one split element
                split_elem(count1)     = iel;                
                tip_enr_pos = find(enrich_node(sctrIGA)==2);
                tip_enr     = sctrIGA(tip_enr_pos);
                
                if isempty(tip_enr)
                    enrich_node(sctrIGA(tNodes)) = 1;
                else
                    hea_enr = setdiff(sctrIGA(tNodes),tip_enr);
                    enrich_node(hea_enr) = 1;
                end
                crack_node(sctrIGA) = iCr;
            end
        end
    end
end
%}



enrich_node = zeros(noCtrPts,1);      % !!!  check in case of no enrichment
crack_node  = zeros(noCtrPts,1);


split_nodes = find(enrich_node == 1);
tip_nodes   = find(enrich_node == 2);

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(node,elementV,'Q4','b-');
cr = plot(xCr(:,1),xCr(:,2),'r-');
set(cr,'LineWidth',3);
n1 = plot(controlPts(split_nodes,1),controlPts(split_nodes,2),'r*');
n2 = plot(controlPts(tip_nodes,1),controlPts(tip_nodes,2),'rs');
set(n1,'MarkerSize',16,'LineWidth',1.07);
set(n2,'MarkerSize',16,'LineWidth',1.07);
axis off  
set(gcf, 'color', 'white');

plot(controlPtsX, controlPtsY,'ro',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',9,'LineWidth',1.0);
            
            
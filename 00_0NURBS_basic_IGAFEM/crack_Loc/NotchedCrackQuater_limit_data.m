% data for edge crack problem with C2 elements

deg = 1;      % degree
cont=deg-1;   % continuity derivative
Rep=deg-cont; % repeated knot

% ref = 7;
Numx = 30;  
Numy = 30;   % luon la so le de vet nut cat qua phtu!!!!!

D = 1 ;      % x dir
L = 1 ;      % y dir
a = 2*L/4;                     % crack length

% coase mesh
[CP,U,V,p,q]=square_coasemesh(D,L,deg);
% p=1;
% q=1;
% U=[0 0 1 1];%m=3 
% V=[0 0 0.5  1 1];%n=3
% %Control Point coordinates
% CP(:,:,1)=[0 0 0  ; L L L];
% CP(:,:,2)=[0 L-a L; 0 L-a L];
% CP(:,:,3)=[1 1 1  ; 1 1 1];
% CP(:,:,4)=[1 1 1  ; 1 1 1];
% [CP,U,V,p,q] = degree_elevate_surf_repeated(p,q,U,V,CP,deg-p,deg-q);

% figure
% plotNURBS_surf_El_CP(p,q,U,V,CP); hold on
% plot_ctrlnet(CP,'ro');view(2)

% [CP,U,V,p,q] = degree_elevate_surf(p,q,U,V,CP,deg-p,deg-q);
% R1 = refinement_vec(U,2*ref);
% R2 = refinement_vec(V,ref);

R1 = refinement_vec_repeated(U,Numx,p-cont);
R2 = refinement_vec_repeated(V,Numy,q-cont);

[CP,uKnot,vKnot] = knot_refine_surf(p,q,U,V,CP,R1,R2);
%b_net=CP(:,:,4);                                        %% chu y chu y

% figure
% plotNURBS_surf_El_CP(p,q,uKnot,vKnot,CP); hold on
% plot_ctrlnet(CP,'ro');view(2)

noPtsX = size(CP,1);       Numx = noPtsX-(p);
noPtsY = size(CP,2);       Numy = noPtsY-(q);

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

weights = reshape(CP(:,:,4),noCtrPts,1);

gcoord(:,1)=reshape(CP(:,:,1),noCtrPts,1);
gcoord(:,2)=reshape(CP(:,:,2),noCtrPts,1);
controlPts=gcoord;

% generate element connectivity ...
% [ien,inn]=genIEN_INN_2D_repeated_Ck(p,q,noPtsX,noPtsY);
[ien,inn]=genIEN_INN_2D(p,q,noPtsX,noPtsY);
element=sort(ien,2);

% generate connectivities and index for physis element
buildVisualizationMesh;
% buildVisualizationMesh_Notched;
% crack data

noCracks = 1;                 % number of cracks

xCrack   = zeros(noCracks,2,2);
xTips    = zeros(noCracks,2);


xCr   = [D 0; D L-a];

xTip  = [D L-a];

seg   = xCr(2,:) - xCr(1,:);   % tip segment
alpha = atan2(seg(2),seg(1));  % inclination angle
QT    =[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];


xCrack(1,:,:) = xCr;
xTips         = xTip;    


% level set computation

numnode   = size(node,1);
numelem   = Numx*Numy;

levelSetCracks

% Choose enriched nodes...
% chooseEnrichedNodes

%
enrich_node = zeros(noCtrPts,1);
crack_node  = zeros(noCtrPts,1); % which crack to which the node enriched

count1 = 0;
count2 = 0;
tol = 1e-8;

tNodes = (p+1)*(p+1);
for i=1:p
    tNodes = [tNodes (p+1)*(p+1)-i*(p+1)]; 
end

% loop over elements
for iel = 1 : numelem
    sctr    = elementV(iel,:);
    sctrIGA = element(iel,:);
    Gnode   = node(sctr,:);
    if and((Gnode(1,1)-Gnode(2,1)~=0) , (Gnode(1,2)-Gnode(4,2)~=0) )
    % loop over cracks
    for iCr = 1 : noCracks        
        phi  = levelSets(iCr,sctr,1);
        psi  = levelSets(iCr,sctr,2);
        
        if ( max(phi)*min(phi) <= tol ) 
            if max(psi) < 0
                count1                 = count1 + 1 ; % ah, one split element
                split_elem(count1)     = iel;                
                tip_enr_pos = find(enrich_node(sctrIGA)==2);
                tip_enr     = sctrIGA(tip_enr_pos);
                
                if isempty(tip_enr)
                    enrich_node(sctrIGA(tNodes)) = 1;
                else
                    hea_enr = setdiff(sctrIGA,tip_enr);
                    enrich_node(hea_enr) = 1;
                end
                crack_node(sctrIGA) = iCr;
            end
        end
        
        if ( max(phi)*min(phi) <= tol ) && (max(psi)*min(psi) < tol)
            count2               = count2 + 1 ; % ah, one tip element
            tip_elem(count2)     = iel;
            enrich_node(sctrIGA) = 2;
            crack_node(sctrIGA)  = iCr;
        end
    end
    end
end
%}
% enrich_node = zeros(noCtrPts,1);      % !!!  check in case of no enrichment
% crack_node  = zeros(noCtrPts,1);

split_nodes = find(enrich_node == 1);
tip_nodes   = find(enrich_node == 2);

% Plot mesh and enriched nodes to check

figure
hold on
plot_mesh(node,elementV,'Q4','b-');
plot(xCr(:,1),xCr(:,2),'r-','LineWidth',3);

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
            
            
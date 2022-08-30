% data for edge crack problem with C0 elements

deg=2; ref=20;  %%% refinement
cont = 0;      % C^0 continuity
% cylinder full
%
p=2;
q=2;
U = [0 0 0 1/4 1/4 2/4 2/4 3/4 3/4 1 1 1]; %m=9
V = [0 0 0 1 1 1]; %n=2
%weights= [ 1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1]';
R_out = 60;
r_in  = 50;
CP(:,:,1)=[R_out (R_out+r_in)/2 r_in;   R_out (R_out+r_in)/2 r_in;   0 0 0;    -R_out -(R_out+r_in)/2 -r_in;...
        -R_out -(R_out+r_in)/2 -r_in;  -R_out -(R_out+r_in)/2 -r_in; 0 0 0;  R_out (R_out+r_in)/2 r_in;  R_out (R_out+r_in)/2 r_in];
     
CP(:,:,2)=[0 0 0;  R_out (R_out+r_in)/2 r_in;   R_out (R_out+r_in)/2 r_in;    R_out (R_out+r_in)/2 r_in;  0 0 0;...
         -R_out -(R_out+r_in)/2 -r_in;  -R_out -(R_out+r_in)/2 -r_in;   -R_out -(R_out+r_in)/2 -r_in;  0 0 0];
CP(:,:,3)=[1 1 1;1 1 1;1 1 1;1 1 1 ;1 1 1;1 1 1;1 1 1;1 1 1 ;1 1 1]*5;
%CP(:,:,4)=[ 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2);1 1];
CP(:,:,4)=[ 1 1 1; 1/sqrt(2) 1/sqrt(2) 1/sqrt(2);1 1 1;1/sqrt(2) 1/sqrt(2) 1/sqrt(2);1 1 1;1/sqrt(2) 1/sqrt(2) 1/sqrt(2);1 1 1;1/sqrt(2) 1/sqrt(2) 1/sqrt(2);1 1 1];
%}

% REFINE uniform=====================================================================
%{
[CP,U,V,p,q] = degree_elevate_surf(p,q,U,V,CP,deg-p,deg-q);
R1 = refinement_vec(U,ref);
R2 = refinement_vec(V,ref);
[CP,U,V] = knot_refine_surf(p,q,U,V,CP,R1,R2);
%}

% REFINE repeated=====================================================================
[CP,U,V,p,q] = degree_elevate_surf_repeated(p,q,U,V,CP,deg-p,deg-q);

R1 = refinement_vec_repeated_p2(U,ref);

R2 = refinement_vec_repeated_p2(V,6);

[CP,uKnot,vKnot] = knot_refine_surf_repeated(p,q,U,V,CP,R1,R2);
%}

plotNURBS_surf_El_CP(p,q,uKnot,vKnot,CP); hold on
plot_ctrlnet(CP,'ro');view(2)

Numx = length(unique(uKnot))-1;
Numy = length(unique(vKnot))-1;

noPtsX = size(CP,1);
noPtsY = size(CP,2);
noCtrPts   = noPtsX * noPtsY;
noDofs     = noCtrPts * 2;

weights    =reshape(CP(:,:,4),noPtsX * noPtsY,1);
gcoord(:,1)=reshape(CP(:,:,1),noPtsX * noPtsY,1);
gcoord(:,2)=reshape(CP(:,:,2),noPtsX * noPtsY,1);
controlPts=gcoord;  %%%%% SAVE controlPts

% generate connectivities and index for control point   
[ien_s,inn]=genIEN_INN_2D_repeated_Ck(p,q,noPtsX,noPtsY);
element=sort(ien_s,2);
ien=ien_s;            %%%% SAVE   ien
% 
% %%%========== remove the overlape control points =================
CP_ovl=[1:noPtsX:noPtsX * noPtsY;
   noPtsX:noPtsX:noPtsX * noPtsY];
% 
% for i=1:size(element,1)
%     index=ien(i,:);
%     for j=1:length(index)
%         [is,loc]=ismember(index(j), CP_ovl(2,:));
%         if is==1
%             ien_s(i,j)=CP_ovl(1,loc);
%         end
%     end
% end
% % element=sort(ien_s,2);
% clear index
% 
% for  j=size(gcoord,1):-1:1
%     is = ismember(j,CP_ovl(2,:));
%     if is ==1
%         gcoord(j,:)=[];
%     end
% end
% controlPts_s=gcoord;

% noCtrPts   = size(controlPts_s,1);
% noDofs      = noCtrPts * 2;

% noPtsX = noPtsX-1; %%%%%%%%%%%%%%%%%%%%%%%% CHU Y
% generate connectivities and index for physis element
buildVisualizationMesh_fullCylinder;

% crack data

noCracks = 1;                 % number of cracks

xCrack   = zeros(noCracks,2,2);
xTips    = zeros(noCracks,2);

a = 0.7*(R_out-r_in);                     % crack length

xCr   = [-r_in 0; -(r_in+a) 0];
xTip  = [-(r_in+a) 0];
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

% %
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
                
        if (min(phi) == 0)  % vet nut di qua phan tu
            %( max(phi)*min(phi) < tol ) % vet nut cat qua phan tu            
            if (max(psi) < 0) && (max(node(sctr,1)) < 0)    %%% crack in left
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


% enrich_node = zeros(noCtrPts,1);      % !!!  check in case of no enrichment
% crack_node  = zeros(noCtrPts,1);


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
            
            
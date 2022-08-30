% data for edge crack problem with C2 elements

p = 3;
q = p;
cont=p-1;   % continuity derivative
Rep=p-cont; % repeated knot

Numx = 19;  
Numy = 19;   % luon la so le de vet nut cat qua phtu!!!!!

D = 1 ;      % x dir
L = 1 ;      % y dir

[CP,U,V,p,q]=square_coasemesh(D,L,p);
R1 = refinement_vec_repeated(U,Numx,p-cont);
R2 = refinement_vec_repeated(V,Numy,q-cont);

[CP,uKnot,vKnot] = knot_refine_surf(p,q,U,V,CP,R1,R2);
%b_net=CP(:,:,4);                                        %% chu y chu y

% % figure
% plotNURBS_surf_El_CP(p,q,uKnot,vKnot,CP); hold on
% plot_ctrlnet(CP,'ro');view(2)

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

noCracks = 2;                 % number of cracks

xCrack   = zeros(noCracks,2,2);
xTips    = zeros(noCracks,2);

a = 2*L/4;                     % crack length

xCr1   = [D/2 0; D/2 (L-a)/2];
xCr2   = [D/2 L; D/2 (L+a)/2];

xTip  = [D/2 (L-a)/2; D/2 (L+a)/2];

seg1   = xCr1(2,:) - xCr1(1,:);   % tip segment
seg2   = xCr2(2,:) - xCr2(1,:);   % tip segment
alpha1 = atan2(seg1(2),seg1(1));  % inclination angle
alpha2 = atan2(seg2(2),seg2(1));  % inclination angle
% QT    =[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];


xCrack(1,:,:) = xCr1;
xCrack(2,:,:) = xCr2;
xTips         = xTip;    


% level set computation

numnode   = size(node,1);
numelem   = size(elementV,1);

levelSetCracks

% Choose enriched nodes...

chooseEnrichedNodes

% enrich_node = zeros(noCtrPts,1);      % !!!  check in case of no enrichment
% crack_node  = zeros(noCtrPts,1);


split_nodes = find(enrich_node == 1);
tip_nodes   = find(enrich_node == 2);

% Plot mesh and enriched nodes to check

figure
hold on
plot_mesh(node,elementV,'Q4','b-');
plot(xCr1(:,1),xCr1(:,2),'r-','LineWidth',3);
plot(xCr2(:,1),xCr2(:,2),'r-','LineWidth',3);

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
            
            
% data for edge crack problem with C2 elements

deg = 3;     ref  = 4;
cont=deg-1;   % continuity derivative
Rep=deg-cont; % repeated knot

p = 2;
q = 2;

U = [0 0 0 1/5 2/5 3/5 4/5 1 1 1 ];    % m=7
V = [0 0 0 1/4 2/4 3/4 1 1 1 ];        % n=6

% load test2.txt;
% gcoordv=test2;
% 
% CP(:,:,1) = reshape(gcoordv(:,1),7,6);
% CP(:,:,2) = reshape(gcoordv(:,2),7,6);
% CP(:,:,3) = reshape(gcoordv(:,3),7,6);
% CP(:,:,4) = reshape(gcoordv(:,4),7,6);


CP(:,:,1) = [      0   0.408241223862591   1.000000000000000   1.000000000000000   0.408241223862590                   0;
                   0   1.000000000000000   1.215995574086253   1.286791611982494   1.000000000000000                   0;
                   0   1.376124652049931   1.498594985511632   1.535644513257927   1.444269822182597                   0;
   2.000000000000000   2.000000000000000   2.000000000000000   2.000000000000000   2.000000000000000   2.000000000000000;
   4.000000000000000   2.584324075253093   2.484725752051881   2.484725752051881   2.548910433004697   4.000000000000000;
   4.000000000000000   3.000000000000000   2.754869851345956   2.674524486090032   3.000000000000000   4.000000000000000;
   4.000000000000000   3.591758776137409   3.000000000000000   3.000000000000000   3.591758776137409   4.000000000000000];

CP(:,:,2) = [3 3.000000000000000 3.605684620961727 4.394315379038273 5.000000000000000 5;
             0 2.000000000000000 3.283908366553260 4.716091633446740 6.000000000000000 8;
             0 1.832355519942581 3.119619690972819 4.880380309027181 6.167644480057419 8;
             0 1.764109279328205 3.000000000000000 5.000000000000000 6.235890720671796 8;
             0 1.832355519942581 3.119619690972819 4.863392018277661 6.167644480057419 8;
             0 2.000000000000000 3.283908366553260 4.716091633446740 6.000000000000000 8;
             3 2.989804580429178 3.605684620961727 4.394315379038273 5.010195419570824 5];
         
CP(:,:,3) = ones(7,6);

CP(:,:,4) = [1.000000  0.844939000000000   0.896626000000000   0.896626000000000   0.844939000000000   1.000000000000000;
   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000;
   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000;
   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000;
   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000;
   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000;
   1.000000000000000   0.844939000000000   0.896626000000000   0.896626000000000   0.844939000000000   1.000000000000000];


% figure      % coasen mesh
% plotNURBS_surf_El_CP(p,q,U,V,CP); hold on
% plot_ctrlnet(CP,'ro');view(2)

[CP,U,V,p,q] = degree_elevate_surf(p,q,U,V,CP,deg-p,deg-q);
R1 = refinement_vec(U,ref);
R2 = refinement_vec(V,ref);

[CP,uKnot,vKnot] = knot_refine_surf(p,q,U,V,CP,R1,R2);

figure
plotNURBS_surf_El_CP(p,q,uKnot,vKnot,CP); hold on
plot_ctrlnet(CP,'ro');view(2)


noPtsX = size(CP,1);      Numx = length(unique(uKnot))-1;
noPtsY = size(CP,2);      Numy = length(unique(vKnot))-1;
noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

weights = reshape(CP(:,:,4),noCtrPts,1);

gcoord(:,1)=reshape(CP(:,:,1),noCtrPts,1);
gcoord(:,2)=reshape(CP(:,:,2),noCtrPts,1);
controlPts=gcoord;

nodes=connecElementQ4(noPtsX-1,noPtsY-1);
%{
 for iel=1:size(nodes,1)
    nd=nodes(iel,:);
    x=gcoord(nd,1);
    y=gcoord(nd,2);    
%    xx=mean(x); 
%    yy=mean(y); 
%      patch(x,y,'w')
     text(x(1),y(1),num2str(nd(1)));
     text(x(2),y(2),num2str(nd(2)));
     text(x(3),y(3),num2str(nd(3)));
     text(x(4),y(4),num2str(nd(4)));
end
%}

% generate connectivities and index for control point   
[ien_temp,inn]= genIEN_INN_2D_repeated_Ck(p,q,noPtsX,noPtsY);

% Loop over elements (knot spans)
tol=1e-8;
ien = [];
%
for e=1:size(ien_temp,1)
    ni = inn(ien_temp(e,1),1);% get NURBS coordinates
    nj = inn(ien_temp(e,1),2);
    if (abs(uKnot(ni)-uKnot(ni+1))>tol)&&(abs(vKnot(nj)-vKnot(nj+1))>tol)
        ien = [ien; ien_temp(e,:)];
    end
end            

element  = sort(ien,2);



% generate connectivities and index for physis element
buildVisualizationMesh;

% crack data

noCracks = 1;                 % number of cracks

% xCrack   = zeros(noCracks,2,2);
% xTips    = zeros(noCracks,2);

a = 0.4;                     % crack length

xCr1   = [3 4; 3-a 4];
% xCr2   = [D/2 L; D/2 (L+a)/2];
% xTip  = [D/2 (L-a)/2; D/2 (L+a)/2];

xTip = [3-a, 4];
seg1   = xCr1(2,:) - xCr1(1,:);   % tip segment
% seg2   = xCr2(2,:) - xCr2(1,:);   % tip segment
alpha1 = atan2(seg1(2),seg1(1));  % inclination angle
% alpha2 = atan2(seg2(2),seg2(1));  % inclination angle
% QT    =[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];


xCrack(1,:,:) = xCr1;
% xCrack(2,:,:) = xCr2;
xTips         = xTip;    


% level set computation

numnode   = size(node,1);
numelem   = size(elementV,1);

levelSetCracks

% Choose enriched nodes...
chooseEnrichedNodes

enrich_node = zeros(noCtrPts,1);      % !!!  check in case of no enrichment
crack_node  = zeros(noCtrPts,1);


split_nodes = find(enrich_node == 1);
tip_nodes   = find(enrich_node == 2);

% Plot mesh and enriched nodes to check

figure
hold on
plot_mesh(node,elementV,'Q4','b-');
plot(xCr1(:,1),xCr1(:,2),'r-','LineWidth',3);
% plot(xCr2(:,1),xCr2(:,2),'r-','LineWidth',3);

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
            
            
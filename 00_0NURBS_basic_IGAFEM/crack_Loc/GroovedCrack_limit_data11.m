% data for edge crack problem with C2 elements

deg = 3;     ref  = 5;
cont=deg-1;   % continuity derivative
Rep=deg-cont; % repeated knot

p = 2;
q = 2;


U = [0 0 0 1/4 2/4 3/4 1 1 1 ];        % n=6
V= [0 0 0 1/5 2/5 3/5 4/5  1 1 1 ];    % m=7

% load test.txt;
% gcoordv=test;
% 
% CP(:,:,1) = reshape(gcoordv(:,1),6,7);
% CP(:,:,2) = reshape(gcoordv(:,2),6,7);
% CP(:,:,3) = reshape(gcoordv(:,3),6,7);
% CP(:,:,4) = reshape(gcoordv(:,4),6,7);


CP(:,:,1) =   [0   0.412467262758882   0.887629714333573   1.033546617199166   0.887629714333573   0.412467262758882 0;
               0   1.000000000000000   1.000000000000000   1.236894146516908   1.000000000000000   1.000000000000000 0;
               0   1.438757175435665   1.507738483519903   1.528921846805270   1.507738483519903   1.438757175435665 0; 
               4   2.561242824564335   2.492261516480097   2.471078153194730   2.492261516480097   2.561242824564335 4;
               4   3.000000000000000   3.000000000000000   2.763105853483092   3.000000000000000   3.000000000000000 4;
               4   3.587532737241118   3.112370285666427   2.966453382800834   3.112370285666427   3.587532737241118 4];


CP(:,:,2) =   [3	3.0246293262	3.4589660969	4	4.5410339031    4.9753706738	5;
               0	1	            2.2460866684	4	5.7539133316	7	            8;
               0	0.8549982887	2.0479755739	4	5.9520244261	7.1450017113	8;
               0	0.8549982887	2.0479755739	4	5.9520244261	7.1450017113	8;
               0	1	            2.2460866684	4	5.7539133316	7	            8;
               3	3.0246293262	3.4589660969	4	4.5410339031	4.9753706738	5];

CP(:,:,3) = ones(6,7);
% 
CP(:,:,4) = [1   0.853553000000000   0.853553000000000   1.000000000000000   0.853553000000000   0.853553000000000 1;
            1   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000 1;
            1   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000 1;
            1   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000 1;
            1   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000 1;
            1   0.853553000000000   0.853553000000000   1.000000000000000   0.853553000000000   0.853553000000000 1];

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
            
            
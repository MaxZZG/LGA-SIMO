% data for center crack problem with C0 elements
% This problem involves two crack tips.

p = 3;      % order
q = p;      % order
cont=0;   % continuity derivative
Rep=p-cont; % repeated knot

w = 1;      % x_dir
D = 2;      % y_dir

Numx=11;  Numy=23;

[CP,U,V,p,q]=square_coasemesh(w,D,p);
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

% generate connectivities and index for physis element
buildVisualizationMesh;

% crack data

a =0.2*w;                     %  crack length
xCr   = [(w-a)/2 D/2; (w+a)/2 D/2];
xTip  = [(w-a)/2 D/2;
         (w+a)/2 D/2];
seg1   = xCr(1,:) - xCr(2,:);   % tip segment 1
seg2   = xCr(2,:) - xCr(1,:);   % tip segment 2

alpha1 = atan2(seg1(2),seg1(1));  % inclination angle
alpha2 = atan2(seg2(2),seg2(1));  % inclination angle

%QT    =[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];


% level set computation

x0  = xCr(1,1); y0 = xCr(1,2);
x1  = xCr(2,1); y1 = xCr(2,2);
t1  = 1/norm(seg1)*seg1;
t2  = 1/norm(seg2)*seg2;

numnode = size(node,1);
numelem = size(elementV,1);

for i = 1 : numnode
    x = node(i,1);
    y = node(i,2);
    l   = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
    phi = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
    ls(i,1) = phi/l;            % normal LS
    ls(i,2) = ([x y]-xTip(1,:))*t1';  % tangent LS tip 1
    ls(i,3) = ([x y]-xTip(2,:))*t2';  % tangent LS tip 1
end

% Choose enriched nodes...

% for one element, if max(phi)*min(phi) < 0
% and max(psi) < 0, then it is a split element
% If max(phi)*min(phi) < 0 and max(psi)*min(psi) < 0, it is
% tip element

% Data structures for elements cut by crack
% Array split_elem contains the number of elements which are completely
% cut by crack. Similarly, array tip_elem stores number of tip element

enrich_node = zeros(noCtrPts,1);
tip_node    = zeros(noCtrPts,1);

crack_node  = ones(noCtrPts,1); % which crack to which the node enriched

count1 = 0;
count2 = 0;

% first Heaviside enriched nodes

for iel = 1 : numelem
    sctr    = elementV(iel,:);
    sctrIGA = element(iel,:);
    phi  = ls(sctr,1);
    psi1 = ls(sctr,2);
    psi2 = ls(sctr,3);
    if ( max(phi)*min(phi) < 0 ) % all elements cut by extended crack
        if max(psi1) < 0 & max(psi2) < 0
            count1 = count1 + 1 ; % ah, one split element
            split_elem(count1)   = iel;
            enrich_node(sctrIGA) = 1;           
        end
    end
end

% then tip enriched nodes, otherwise, some tip enriched
% nodes can be overwritten by H enriched ones
tol=1e-10;
for iel = 1 : numelem
    sctr    = elementV(iel,:);
    sctrIGA = element(iel,:);
    phi  = ls(sctr,1);
    psi1 = ls(sctr,2);
    psi2 = ls(sctr,3);
    if ( max(phi)*min(phi) < tol ) % all elements cut by extended crack
        if max(psi1)*min(psi1) < tol
            count2 = count2 + 1 ; % ah, one tip 1 element
            tip_elem(count2)      = iel;
            enrich_node(sctrIGA)  = 2;
            tip_node(sctrIGA)     = 1;
        elseif max(psi2)*min(psi2) < tol
            count2 = count2 + 1 ; % ah, one tip 2 element
            tip_elem(count2)     = iel;
            enrich_node(sctrIGA) = 2;            
            tip_node(sctrIGA)    = 2;
        end
    end
end

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
set(n1,'MarkerSize',14,'LineWidth',1.07);
set(n2,'MarkerSize',14,'LineWidth',1.07);
axis off  
set(gcf, 'color', 'white');

plot(controlPtsX, controlPtsY,'ro',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',8,'LineWidth',1.0);
            
            
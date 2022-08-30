%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended Isogeometric analysis for two dimensional linear elastic
% interfacial fracture mechanics problems.
%
% Interfacial center crack plate in tension (half model)
% Level sets are used to detect enriched nodes
% They are not used in computing Heaviside and branch functions.
%
% Vinh Phu Nguyen, April 2012, Hue Vietnam
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/
addpath ../nurbs-geopdes/inst/
addpath ../nurbs-util/

close all
clear variables

global p q controlPts weights element xCrack xTips crack_node ...
       ep CHI jDomainFac weakEnrFunc


nu1          = 0.2571;  % Poisson ratio
E2           = 100;  % of mat2
E1           = 22*E2;  % Young modulus of mat1
nu2          = 0.3;
stressState = 'PLANE_STRAIN';

weakEnrFunc     = 1; %1=Moes func 2=absolute signed dis
layerEnrichment = 0; %0(off)
tipEnrId        = 4; %=4:13 enr. functions (12 bi-mat branch + Moes func)
                     %=5:12 enr. functions (12 bi-mat branch )
                     %=6:5 enr. functions (4 homo. branch + Moes func)
                     %=2:4 homogeneous branch functions
inclusionEnrId  = 3; % =3: weak discontinuity enrichment
                     % =0: no enrichment
jDomainFac      = 2; % J-integral domain = jDomainFac*area tip elem

vtuFile     = '../results/edgeBiMatCrackXIGA';

% Elasticity matrices

Cm = elasticityMatrix(E1,nu1,stressState);
Ci = elasticityMatrix(E2,nu2,stressState);

% Constant in the 12 tip enrichment functions

G1 = E1/2/(1+nu1);       % Shear modulus for mat1
G2 = E2/2/(1+nu2);       % Shear modulus for mat2

% Kosolov constants

if strcmp(stressState,'PLANE_STRESS')
    k1 = (3-nu1)/(1+nu1);
    k2 = (3-nu2)/(1+nu2);
else
    k1 = 3-4*nu1;
    k2 = 3-4*nu2;
end

b   = (G1*(k2-1)-G2*(k1-1))/(G1*(k2+1)+G2*(k1+1));    % Second Dundur's parameter
ep  = 1/(2*pi)*log((1-b)/(1+b));                      % Material constant

% cos(-ep*log(2*a)) - 2*ep*sin(-ep*log(2*a))
% sin(-ep*log(2*a)) + 2*ep*cos(-ep*log(2*a))

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lineInterfaceCrackCenterC2Data 


% find boundary nodes for boundary conditions

leftEdge     =  find(controlPts(:,1)==0); %left edge
rightEdge    =  find(controlPts(:,1)==W); %right edge
topNodes     =  find(controlPts(:,2)==2*W); %top edge
botNodes     =  find(controlPts(:,2)==0); %bottom edge

fixedXNodes  =  [leftEdge; rightEdge];
fixedYNodes = botNodes;

% build boundary mesh for force vector computation

topEdgeMesh    = zeros(noElemsU,p+1);
botEdgeMesh    = zeros(noElemsU,p+1);

for i=1:noElemsU
    topEdgeMesh(i,:) = topNodes(i:i+p);
    botEdgeMesh(i,:) = botNodes(i:i+p);
end

% Each split node is enriched by ONE function, H(x)
% Each tip node is enriched by FOUR functions, B_i(x)
% then, the total dofs is :
% total dof = numnode*nsd + numsplitnode*1*nsd + numstipnode*4*nsd
% here, two dimension, nsd = 2

% Each split node is enriched by ONE function, H(x)
% Each tip node is enriched by FOUR functions, B_i(x)
% then, the total dofs is :
% total dof = numnode*nsd + numsplitnode*1*nsd + numstipnode*4*nsd
% here, two dimension, nsd = 2

noDofs = noDofs + size(split_nodes,1)*1*2 + ...
                  size(tip_nodes,1)*4*2 + size(inc_nodes,1)*1*2 +...
                  size(itip_nodes,1)*13*2 + ... % 12 bi-mat branch fns+Moes fnc
                  size(iTIp_nodes,1)*5*2 + ...  % 4 homo. branch fns+Moes fnc
                  size(iTip_nodes,1)*12*2;      % % 12 bi-mat branch fns

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

% We use fictitious nodes/control points to handle these
% additional dofs. At a H(x) enriched node, we add one fantom node and
% at tip enriched node, four fantom nodes are added. These fictitious nodes
% are numbered from the total number of true nodes, ie, from numnode+1 ...

[pos] = computePosOfEnrichedNodes(noCtrPts,enrich_node);


% essential boundary conditions

uFixed     = zeros(size(fixedXNodes))';
vFixed     = zeros(size(fixedYNodes))';

udofs      = 2*fixedXNodes-1;
vdofs      = 2*fixedYNodes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assembling system of equation
% Stiffness matrix and external force vector

noGPs         = p+1; % # of Gauss points along one direction
noGP1         = p+1;

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

gps = [];

% Loop over elements (knot spans)

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);          %  element scatter vector
    sctrV  = elementV(e,:);         %  element scatter vector Q4
    levelS = CHI(1,sctr);     
    nn     = length(sctr);
    pts    = controlPts(sctr,:);
    
    % -----------------------------------------------
    % Choose Gauss quadrature rules for elements
    
    [W,Q] = gaussForEnrichedElement(e,noGPs,sctrV,levelSets,xTip,split_elem, splitElems,tip_elem,...
                                   itip_nodes,iTip_nodes);
    
    % Determine the position in the global matrix K
    
    sctrB = assembly(e,enrich_node,pos);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % compute derivative of basis functions w.r.t parameter coord
        
        [N dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');
        
        % B matrix
        
        [B,J1] = BMatrixXIGA(Xi,Eta,e,enrich_node,N,dRdxi,dRdeta);
        
        % Stiffness matrix
        
        levelset = dot(N,levelS);
        
        if levelset >= 0
            C = Cm;
        else
            C = Ci;
        end
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B'*C*B*W(gp)*J1*J2;
        
        x         = N * pts; % global coord of GP
        gps       = [gps;x];
    end
end

% Computing external force

[W1,Q1] = quadrature(noGP1, 'GAUSS', 1 );

% Loop over elements along top edge

bndPoints  = controlPts(topNodes,:);

sigma0 = 1;

for e=1:noElemsU
    xiE   = elRangeU(e,:); % [xi_i,xi_i+1]
    conn  = elConnU(e,:);
    pts   = bndPoints(conn,:);
    
    sctrx = topEdgeMesh(e,:);
    sctry = sctrx*2;
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N dNdxi] = NURBS1DBasisDers(Xi,p,uKnot,weights);
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob1   = dNdxi * pts;
        J1       = norm (jacob1);
        
        f(sctry) = f(sctry) + N' * sigma0 * J1 * J2 * wt;
    end
end

bndPoints  = controlPts(botNodes,:);

% for e=1:noElemsU
%     xiE   = elRangeU(e,:); % [xi_i,xi_i+1]
%     conn  = elConnU(e,:);
%     pts   = bndPoints(conn,:);
%     
%     sctrx = botEdgeMesh(e,:);
%     sctry = sctrx*2;
%     
%     % loop over Gauss points
%     for gp=1:size(W1,1)
%         xi      = Q1(gp,:);
%         wt      = W1(gp);
%         Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
%         J2      = 0.5 * ( xiE(2) - xiE(1) );
%         
%         [N dNdxi] = NURBS1DBasisDers(Xi,p,uKnot,weights);
%         
%         % compute the jacobian of physical and parameter domain mapping
%         % then the derivative w.r.t spatial physical coordinates
%         
%         jacob1   = dNdxi * pts;
%         J1       = norm (jacob1);
%         
%         f(sctry) = f(sctry) - N' * sigma0 * J1 * J2 * wt;
%     end
% end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

applyBC

% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux = U(1:2:2*noCtrPts);
Uy = U(2:2:2*noCtrPts);

% compute stresses and displacements at mesh vertices
disp([num2str(toc),'  COMPUTING STRESSES '])

plot(gps(:,1),gps(:,2),'+');

plotStressXIGAMultiMats

% plot stress

% stressComp=2;
% figure
% clf
% plot_field(node,elementV,'Q4',stress(:,:,stressComp));
% hold on
% colorbar
% title('Stress in x direction')
% axis off
% %plot_mesh(node,elementV,'Q4','g.-');
% 
% figure
% clf
% fac=0.1;
% plot_field(node+fac*[dispX dispY],elementV,'Q4',dispY);
% hold on
% colorbar
% title('Displacement in y direction')
% axis off

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear disp;

disp([num2str(toc),'  SIFs COMPUTATION'])

computeInteractionIntegralsBiMatCracks

F1 = KI/sigma0/sqrt(pi*a)
F2 = KII/sigma0/sqrt(pi*a)

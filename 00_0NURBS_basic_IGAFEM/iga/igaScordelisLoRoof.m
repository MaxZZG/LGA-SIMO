%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for three dimensional elasticity problems.
%
% Scordelis-Lo roof under gravity loading
% 1/4 model is analysed.
%
% Vinh Phu Nguyen,
% Johns Hopkins University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../nurbs-util/
addpath ../integration/

close all
clear variables

global p q r uKnot vKnot wKnot

E0           = 4.32e8;  % Young's modulus
nu0          = 0.0;  % Poisson’s ratio


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE COMPLIANCE MATRIX
C=zeros(6,6);
C(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
                                  nu0 1-nu0 nu0;
                                  nu0 nu0 1-nu0];
C(4:6,4:6)=E0/(1+nu0)*eye(3);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
scordelisLoCkData % new data using NURBS toolbox
                      % that allows k-refinement 

plotMesh3(controlPts,weights, uKnot,vKnot,wKnot,...
                   p,q,r,40,'r-','try.eps');
               
noGPs      = 4; % # of Gauss points along one direction

noCtrPts   = noPtsX * noPtsY * noPtsZ;
noDofs     = noCtrPts * 3;
g          = -90/t;

%% Boundary nodes
% find boundary nodes for boundary conditions

eps = 1e-14;

rigidNodes  = find(controlPts(:,3)==0)'; 

zConsNodes  = find(abs(controlPts(:,3)-L)<=eps)';
xConsNodes  = find(controlPts(:,1)==0)';

xConsNodes  = [xConsNodes rigidNodes];
yConsNodes  = rigidNodes;

xConsNodes  = unique(xConsNodes);
zConsNodes  = unique(zConsNodes);

% essential boundary conditions

uFixed     = zeros(size(xConsNodes));
vFixed     = zeros(size(yConsNodes));
wFixed     = zeros(size(zConsNodes));

udofs      = xConsNodes;             % global indecies of the fixed x disps
vdofs      = yConsNodes+noCtrPts;    % global indecies of the fixed y disps
wdofs      = zConsNodes+2*noCtrPts;  % global indecies of the fixed z disps


%% build connectivity ...

generateIGA3DMesh

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
%[W,Q]=quadrature(  noGPs, 'GAUSS', 3 ); 
[W,Q]=gaussianQuadNURBS(p+1,q+1,r+1); 

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over elements (knot spans)

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    idw    = index(e,3);
    
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    zetaE  = elRangeW(idw,:); % [zeta_k,zeta_k+1]
 
    sctr   = element(e,:);          %  element scatter vector
    sctrB  = [sctr sctr+noCtrPts sctr+2*noCtrPts]; % scatters a B matrix
    nn     = length(sctr);
    pts    = controlPts(sctr,:);
    B      = zeros(6,3*nn);
    sctry  = sctr + noCtrPts;  
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,  pt(1));
        Eta     = parent2ParametricSpace(etaE, pt(2));
        Zeta    = parent2ParametricSpace(zetaE,pt(3));
        J2      = jacobianPaPaMapping3d (xiE,etaE,zetaE);
        
        % compute derivative of basis functions w.r.t parameter coord
        
        [N, dRdxi, dRdeta, dRdzeta] = NURBS3DBasisDers([Xi;Eta;Zeta],...
                                   p,q,r,uKnot,vKnot,wKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates                
        % Jacobian matrix
             
        jacob      = pts'*[dRdxi' dRdeta' dRdzeta'];
        J1         = det(jacob);                        
        
        dRdx       = [dRdxi' dRdeta' dRdzeta']/jacob;
        
        % B matrix
        
        B          = strainDispMatrix3d(nn,dRdx);
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        f(sctry)       = f(sctry)       + N' * g * J1 * J2 * wt;
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J1 * J2 * wt;
    end
end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

bcwt=mean(diag(K)); % a measure of the average  size of an element in K
                    % used to keep the  conditioning of K matrix

f=f-K(:,udofs)*uFixed';  % modify the  force vector
f=f-K(:,vdofs)*vFixed';
f=f-K(:,wdofs)*wFixed';

f(udofs) = bcwt*uFixed;
f(vdofs) = bcwt*vFixed;
f(wdofs) = bcwt*wFixed;

K(udofs,:)=0;  % zero out the rows and  columns of the K matrix
K(vdofs,:)=0;
K(wdofs,:)=0;

K(:,udofs)=0;
K(:,vdofs)=0;
K(:,wdofs)=0;

K(udofs,udofs)=bcwt*speye(length(udofs));  % put ones*bcwt on the diagonal
K(vdofs,vdofs)=bcwt*speye(length(vdofs));
K(wdofs,wdofs)=bcwt*speye(length(wdofs));

% SOLVE SYSTEM

disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([num2str(toc),'  POST-PROCESSING'])

Ux    = U(1:noCtrPts);
Uy    = U(noCtrPts+1:2*noCtrPts);
Uz    = U(2*noCtrPts+1:noDofs);

num_uy = min(Uy)
exa_uy = -0.3024

vtsFile = '../results/scordelisLoRoof';
plotStress3d






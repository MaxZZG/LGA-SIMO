function phi = StreamLines_ls(X,T,elem,velo,npx,npy,vect,fig_num)
% phi = StreamLines_ls(X,T,elem,velo,npx,npy,vect,fig_num)
% Computation of a function phi such that its level sets
% are streamlines of the velocity field velo
%                   phi_ii = u_y-v_x  
%

if nargin < 6
    error('Too few input arguments')
elseif nargin == 6
    vect = [];
    fig_num = [];
elseif nargin == 7
    fig_num = [];
end

nx = npx-1; 
ny = npy-1; 
npt = npx*npy;

% Matrices
L = CreLaplace(X,T,elem); 
f = CreVect_phi(X,T,elem,velo); 

% Boundary conditions
numnp = size(L,1);
[Accd, bccd] =  BC_streamlines(nx,ny,numnp);
nDir = size(Accd,1);

% Entire matrices (including boundary conditions)
Ltot = [L Accd'; Accd zeros(nDir)];
ftot = [f;bccd];

% Solution
phi = Ltot\ftot;

if isempty(vect)
    min_phi = min(phi);
    max_phi = max(phi);
    aux = max_phi - min_phi;
    vect1 = linspace(min_phi,0,10);
    vect2 = linspace(0,-min_phi,5);
    vect3 = linspace(aux,max_phi,15);
    vect4 = linspace(-min_phi,aux,70);
    vect = [vect1,vect2,vect3,vect4];
end


% Postprocess: streamlines
aux_x = ones(npy,1)*X(1:npx,1)';
aux_y = X(1:npx:npt,2)*ones(1,npx);
aux_phi = reshape(phi(1:npt),npx,npy)';

if isempty(fig_num)
    figure('Name','Streamlines computed as level sets','NumberTitle','off'); clf; 
else
    figure(fig_num); clf;
end
contour(aux_x,aux_y,aux_phi,vect)
box on; axis equal tight;
set(gca,'FontSize',12,...
    'XTick',[0:0.25:1], 'YTick',[0:0.25:1]);
axis([0,1,0,1])

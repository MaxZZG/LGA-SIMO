function [Accd, bccd, velo] =  BC_Cavity(nx,ny,nunk)
% [Accd, bccd, velo] =  BC_Cavity(nx,ny,nunk)
% This function creates matrices Accd and bccd to impose Dirichlet 
% boudary conditions for the Cavity flow problem using Lagrange 
% multipliers method.
% Velocity is prescibed in the whole boundary and vertical component
% is zero everywhere whereas the horizontal one is zero everywhere 
% but the top, where it is equal one. 
% Besides, the function provides a velocity field velo which verifies boundary
% conditions and can be used as an initial guess for a nonlinear analysis.
%
% Input:
%   nx,ny:  number of elements on each direction
%   nunk:   number of degrees of freedom for the velocity field
%

% Nodes on the Dirichlet boundary
nodesDy0 = [1 : nx+1]';                             % y = 0
nodesDx1 = [2*(nx+1) : nx+1 : (ny+1)*nx]';          % x = 1
nodesDx0 = [(ny-1)*(nx+1)+1 : -(nx+1) : nx+2]';     % x = 0
nodesDy1 = [(ny+1)*(nx+1) : -1 : ny*(nx+1)+1]';     % y = 1
nodesD   = [nodesDy0;nodesDx1;nodesDx0;nodesDy1];

% Imposed boundary conditions
Cy0 = [reshape([2*nodesDy0'-1;            2*nodesDy0']              ,2*size(nodesDy0,1),1), ...
       reshape([zeros(1,size(nodesDy0,1));zeros(1,size(nodesDy0,1))],2*size(nodesDy0,1),1)];
Cx1 = [reshape([2*nodesDx1'-1;            2*nodesDx1']              ,2*size(nodesDx1,1),1), ...
       reshape([zeros(1,size(nodesDx1,1));zeros(1,size(nodesDx1,1))],2*size(nodesDx1,1),1)];
Cx0 = [reshape([2*nodesDx0'-1;            2*nodesDx0']              ,2*size(nodesDx0,1),1), ...
       reshape([zeros(1,size(nodesDx0,1));zeros(1,size(nodesDx0,1))],2*size(nodesDx0,1),1)];
Cy1 = [reshape([2*nodesDy1'-1;            2*nodesDy1']              ,2*size(nodesDy1,1),1), ...
      reshape([ones(1,size(nodesDy1,1)) ;zeros(1,size(nodesDy1,1))],2*size(nodesDy1,1),1)];
C = [Cy0; Cx1; Cx0; Cy1];

% Boundary conditions' matrix
nDir = size(C,1);
Accd = zeros(nDir,nunk);
Accd(:,C(:,1)) = eye(nDir);
% Boundary conditions' vector
bccd = C(:,2);

% Initial velocity: null field everywhere but on the boundary, 
% where essential boundary conditions are satisfied.
velo = zeros(nunk/2,2);
velo(nodesD,:) = reshape(C(:,2),2,size(nodesD,1))';


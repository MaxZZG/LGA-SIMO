function [N, dN, ddN, coord, J]=nurbshape3da(e,u,v,w,u_knot,v_knot,w_knot,b_net,p,q,r,~,~,~,element_nod,coord_ijk)

%calculate the shape function and second derivatives 


nsd = 3; 
% ------------------------------------------------------------------
%     get nurbs coordinates forml local node 1;

ni = coord_ijk(element_nod(e,end),1);
nj = coord_ijk(element_nod(e,end),2);
nk = coord_ijk(element_nod(e,end),3);


mcp = 2;
ncp = 2;
pcp = 2;

%     evaluate 1d size functions and derivatives each direction;
% M=derBasisFuns(ni,p,mcp,u,u_knot); % calculate in u direction
% P=derBasisFuns(nj,q,ncp,v,v_knot); % calculate in v direction
% Q=derBasisFuns(nk,r,pcp,w,w_knot); % calculate in w direction

M=dersbasisfuns2(ni,p,mcp,u,u_knot); % calculate in u direction
P=dersbasisfuns2(nj,q,ncp,v,v_knot); % calculate in v direction
Q=dersbasisfuns2(nk,r,pcp,w,w_knot); % calculate in w direction

wgts= zeros((p+1)*(q+1)*(r+1), 1);
cpts = zeros((p+1)*(q+1)*(r+1), nsd);
B = zeros(1, (p+1)*(q+1)*(r+1));
dBxi = zeros(1, (p+1)*(q+1)*(r+1));
dBeta = zeros(1, (p+1)*(q+1)*(r+1));
dBzeta = zeros(1, (p+1)*(q+1)*(r+1));
ddBxi = zeros(1, (p+1)*(q+1)*(r+1));
ddBeta = zeros(1, (p+1)*(q+1)*(r+1));
ddBzeta = zeros(1, (p+1)*(q+1)*(r+1));
ddBxieta = zeros(1, (p+1)*(q+1)*(r+1));
ddBxizeta = zeros(1, (p+1)*(q+1)*(r+1));
ddBetazeta = zeros(1, (p+1)*(q+1)*(r+1));

N = zeros(1, (p+1)*(q+1)*(r+1));
dN = zeros(3, (p+1)*(q+1)*(r+1));
ddN = zeros(6, (p+1)*(q+1)*(r+1));

% Form tensor products
icount = 0;

for k = 1:r+1 
    for j = 1:q+1
        for i = 1:p+1
            icount = icount+1;

            wgts(icount) = b_net(ni+i-p-1,nj+j-q-1,nk+k-r-1,nsd+1);
            cpts(icount, 1) = b_net(ni+i-p-1,nj+j-q-1,nk+k-r-1,1);
            cpts(icount, 2) = b_net(ni+i-p-1,nj+j-q-1,nk+k-r-1,2);
            cpts(icount, 3) = b_net(ni+i-p-1,nj+j-q-1,nk+k-r-1,3);
        %   basis functions;    
            B(icount) = M(1,i)*P(1,j)*Q(1,k);
            dBxi(icount) = M(2,i)*P(1,j)*Q(1,k);
            dBeta(icount) = M(1,i)*P(2,j)*Q(1,k);          
            dBzeta(icount) = M(1,i)*P(1,j)*Q(2,k);
            ddBxi(icount) = M(3,i)*P(1,j)*Q(1,k);
            ddBeta(icount) = M(1,i)*P(3,j)*Q(1,k);
            ddBzeta(icount) = M(1,i)*P(1,j)*Q(3,k);
            ddBxieta(icount) = M(2,i)*P(2,j)*Q(1,k);
            ddBxizeta(icount) = M(2,i)*P(1,j)*Q(2,k);
            ddBetazeta(icount) = M(1,i)*P(2,j)*Q(2,k);
        end
    end 
end

    
% Multiply each B-spline function with corresponding weight
N(1,:) = B .* wgts';
dN(1,:) = dBxi .* wgts';
dN(2,:) = dBeta .* wgts';
dN(3,:) = dBzeta .* wgts';
ddN(1,:) = ddBxi .* wgts';
ddN(2,:) = ddBeta .* wgts';
ddN(3,:) = ddBzeta .* wgts';
ddN(4,:) = ddBxieta .* wgts';
ddN(5,:) = ddBxizeta .* wgts';
ddN(6,:) = ddBetazeta .* wgts';

% Compute the sums of B-spline functions
w_sum = sum(N(1,:));
dw_xi = sum(dN(1,:));
dw_eta = sum(dN(2,:));
dw_zeta = sum(dN(3,:));
d2w_xi = sum(ddN(1,:));
d2w_eta = sum(ddN(2,:));
d2w_zeta = sum(ddN(3,:));
d2w_xieta = sum(ddN(4,:));
d2w_xizeta = sum(ddN(5,:));
d2w_etazeta = sum(ddN(6,:));

% Compute NURBS basis functions and its first and second derivatives in
% local coordinates
ddN(1,:) = ddN(1,:)/w_sum - (2*dN(1,:)*dw_xi + N*d2w_xi)/w_sum^2 + 2*N*dw_xi^2/w_sum^3; %dxidxi derivative
ddN(2,:) = ddN(2,:)/w_sum - (2*dN(2,:)*dw_eta + N*d2w_eta)/w_sum^2 + 2*N*dw_eta^2/w_sum^3; %detadeta derivative
ddN(3,:) = ddN(3,:)/w_sum - (2*dN(3,:)*dw_zeta + N*d2w_zeta)/w_sum^2 + 2*N*dw_zeta^2/w_sum^3; %dzetadzeta derivative
ddN(4,:) = ddN(4,:)/w_sum - (dN(1,:)*dw_eta + dN(2,:)*dw_xi + N*d2w_xieta)/w_sum^2 + ...
    2*N*dw_xi*dw_eta/w_sum^3; %dxideta derivative
ddN(5,:) = ddN(5,:)/w_sum - (dN(1,:)*dw_zeta + dN(3,:)*dw_xi + N*d2w_xizeta)/w_sum^2 + ...
    2*N*dw_xi*dw_zeta/w_sum^3; %dxidzeta derivative
ddN(6,:) = ddN(6,:)/w_sum - (dN(2,:)*dw_zeta + dN(3,:)*dw_eta + N*d2w_etazeta)/w_sum^2 + ...
    2*N*dw_eta*dw_zeta/w_sum^3; %detadzeta derivative

dN(1,:) = dN(1,:)/w_sum - N*dw_xi/w_sum^2;
dN(2,:) = dN(2,:)/w_sum - N*dw_eta/w_sum^2;
dN(3,:) = dN(3,:)/w_sum - N*dw_zeta/w_sum^2;

N = N/w_sum;


% calculate coordinates in physical space
coord = N*cpts;

% Compute Jacobian matrix
dxdxi = dN*cpts;
dxdxi = dxdxi';
J = det(dxdxi);


% Set up the second derivatives matrix and the matrix of squared first derivatives
d2xdxi2 = ddN*cpts;

dxdxi2 = [ dxdxi(1,1)^2            dxdxi(1,2)^2             dxdxi(1,3)^2            dxdxi(1,1)*dxdxi(1,2)                       dxdxi(1,1)*dxdxi(1,3)                       dxdxi(1,2)*dxdxi(1,3)
           dxdxi(2,1)^2            dxdxi(2,2)^2             dxdxi(2,3)^2            dxdxi(2,1)*dxdxi(2,2)                       dxdxi(2,1)*dxdxi(2,3)                       dxdxi(2,2)*dxdxi(2,3)
           dxdxi(3,1)^2            dxdxi(3,2)^2             dxdxi(3,3)^2            dxdxi(3,1)*dxdxi(3,2)                       dxdxi(3,1)*dxdxi(3,3)                       dxdxi(3,2)*dxdxi(3,3)
           2*dxdxi(1,1)*dxdxi(2,1) 2*dxdxi(1,2)*dxdxi(2,2)  2*dxdxi(1,3)*dxdxi(2,3) dxdxi(1,1)*dxdxi(2,2)+dxdxi(1,2)*dxdxi(2,1) dxdxi(1,1)*dxdxi(2,3)+dxdxi(1,3)*dxdxi(2,1) dxdxi(1,2)*dxdxi(2,3)+dxdxi(1,3)*dxdxi(2,2)
           2*dxdxi(1,1)*dxdxi(3,1) 2*dxdxi(1,2)*dxdxi(3,2)  2*dxdxi(1,3)*dxdxi(3,3) dxdxi(1,1)*dxdxi(3,2)+dxdxi(1,2)*dxdxi(3,1) dxdxi(1,1)*dxdxi(3,3)+dxdxi(1,3)*dxdxi(3,1) dxdxi(1,2)*dxdxi(3,3)+dxdxi(1,3)*dxdxi(3,2)
           2*dxdxi(2,1)*dxdxi(3,1) 2*dxdxi(2,2)*dxdxi(3,2)  2*dxdxi(2,3)*dxdxi(3,3) dxdxi(2,1)*dxdxi(3,2)+dxdxi(2,2)*dxdxi(3,1) dxdxi(2,1)*dxdxi(3,3)+dxdxi(2,3)*dxdxi(3,1) dxdxi(2,2)*dxdxi(3,3)+dxdxi(2,3)*dxdxi(3,2) ];

% Solve for first derivatives in global coordinates 

dN = dxdxi'\dN;

%dN = pinv(dxdxi')*dN;

% Solve for second derivatives in global coordinates 
ddN = dxdxi2'\(ddN - d2xdxi2*dN);
%ddN = pinv(dxdxi2')*(ddN - d2xdxi2*dN);


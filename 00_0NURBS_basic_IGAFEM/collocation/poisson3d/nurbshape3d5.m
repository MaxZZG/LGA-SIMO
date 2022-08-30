function [N, dN, ddN]=nurbshape3d5(M, P, Q, p, q, r, wgts)
%calculates the shape functions, first and second derivatives with respect
%to the parameter space

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


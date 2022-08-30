function [N, dN, ddN]=nurbshape2d5(M, P, p, q, wgts)
%calculates the shape functions, first and second derivatives with respect
%to the parameter space

B = zeros(1, (p+1)*(q+1));
dBxi = zeros(1, (p+1)*(q+1));
dBeta = zeros(1, (p+1)*(q+1));
ddBxi = zeros(1, (p+1)*(q+1));
ddBeta = zeros(1, (p+1)*(q+1));
ddBxieta = zeros(1, (p+1)*(q+1));

N = zeros(1, (p+1)*(q+1));
dN = zeros(2, (p+1)*(q+1));
ddN = zeros(3, (p+1)*(q+1));

% Form tensor products
icount = 0;


for j = 1:q+1
    for i = 1:p+1
        icount = icount+1;

    %   basis functions;    
        B(icount) = M(1,i)*P(1,j);
        dBxi(icount) = M(2,i)*P(1,j);
        dBeta(icount) = M(1,i)*P(2,j);        
        ddBxi(icount) = M(3,i)*P(1,j);
        ddBeta(icount) = M(1,i)*P(3,j);        
        ddBxieta(icount) = M(2,i)*P(2,j);        
    end
end 

    
% Multiply each B-spline function with corresponding weight
N(1,:) = B .* wgts';
dN(1,:) = dBxi .* wgts';
dN(2,:) = dBeta .* wgts';
ddN(1,:) = ddBxi .* wgts';
ddN(2,:) = ddBxieta .* wgts';
ddN(3,:) = ddBeta .* wgts';

% Compute the sums of B-spline functions
w_sum = sum(N(1,:));
dw_xi = sum(dN(1,:));
dw_eta = sum(dN(2,:));
d2w_xi = sum(ddN(1,:));
d2w_xieta = sum(ddN(2,:));
d2w_eta = sum(ddN(3,:));

% Compute NURBS basis functions and its first and second derivatives in
% local coordinates
ddN(1,:) = ddN(1,:)/w_sum - (2*dN(1,:)*dw_xi + N*d2w_xi)/w_sum^2 + 2*N*dw_xi^2/w_sum^3; %dxidxi derivative
ddN(3,:) = ddN(3,:)/w_sum - (2*dN(2,:)*dw_eta + N*d2w_eta)/w_sum^2 + 2*N*dw_eta^2/w_sum^3; %detadeta derivative
ddN(2,:) = ddN(2,:)/w_sum - (dN(1,:)*dw_eta + dN(2,:)*dw_xi + N*d2w_xieta)/w_sum^2 + ...
    2*N*dw_xi*dw_eta/w_sum^3; %dxideta derivative

dN(1,:) = dN(1,:)/w_sum - N*dw_xi/w_sum^2;
dN(2,:) = dN(2,:)/w_sum - N*dw_eta/w_sum^2;

N = N/w_sum;


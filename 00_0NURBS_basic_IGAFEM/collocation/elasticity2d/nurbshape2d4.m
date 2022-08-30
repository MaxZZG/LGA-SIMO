function [N, dN]=nurbshape2d4(M, P, p, q, wgts)

%calculate the shape function and derivatives with respect to the parameter
%space

B = zeros(1, (p+1)*(q+1));
dBxi = zeros(1, (p+1)*(q+1));
dBeta = zeros(1, (p+1)*(q+1));

N = zeros(1, (p+1)*(q+1));
dN = zeros(2, (p+1)*(q+1));


% Form tensor products
icount = 0;
for  j= 1:q+1 
    for i = 1:p+1        
        icount = icount+1;
    %   basis functions;    
        B(icount) = M(1,i)*P(1,j);
        dBxi(icount) = M(2,i)*P(1,j);
        dBeta(icount) = M(1,i)*P(2,j);                              
    end 
end

    
% Multiply each B-spline function with corresponding weight
N(1,:) = B .* wgts';
dN(1,:) = dBxi .* wgts';
dN(2,:) = dBeta .* wgts';

% Compute the sums of B-spline functions
w_sum = sum(N(1,:));
dw_xi = sum(dN(1,:));
dw_eta = sum(dN(2,:));

% Compute NURBS basis functions and its first derivatives in
% local coordinates

dN(1,:) = dN(1,:)/w_sum - N*dw_xi/w_sum^2;
dN(2,:) = dN(2,:)/w_sum - N*dw_eta/w_sum^2;

N = N/w_sum;

function [N, dN]=nurbshape3d4(M, P, Q, p, q, r, wgts)

%calculate the shape function and derivatives with respect to the parameter
%space

B = zeros(1, (p+1)*(q+1)*(r+1));
dBxi = zeros(1, (p+1)*(q+1)*(r+1));
dBeta = zeros(1, (p+1)*(q+1)*(r+1));
dBzeta = zeros(1, (p+1)*(q+1)*(r+1));


N = zeros(1, (p+1)*(q+1)*(r+1));
dN = zeros(3, (p+1)*(q+1)*(r+1));


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
        end
    end 
end

% TO CHECK: is the order the same as expected for when wgts are not all 1s?
% nument = (p+1)*(q+1)*(r+1);
% B = reshape(kron(M(1,:)'*P(1,:),Q(1,:)'), 1, nument);
% dBxi = reshape(kron(M(2,:)'*P(1,:),Q(1,:)'), 1, nument);
% dBeta = reshape(kron(M(1,:)'*P(2,:),Q(1,:)'), 1, nument);
% dBzeta = reshape(kron(M(1,:)'*P(1,:),Q(2,:)'), 1, nument);


    
% Multiply each B-spline function with corresponding weight
N(1,:) = B .* wgts';
dN(1,:) = dBxi .* wgts';
dN(2,:) = dBeta .* wgts';
dN(3,:) = dBzeta .* wgts';

% Compute the sums of B-spline functions
w_sum = sum(N(1,:));
dw_xi = sum(dN(1,:));
dw_eta = sum(dN(2,:));
dw_zeta = sum(dN(3,:));


% Compute NURBS basis functions and its first derivatives in
% local coordinates

dN(1,:) = dN(1,:)/w_sum - N*dw_xi/w_sum^2;
dN(2,:) = dN(2,:)/w_sum - N*dw_eta/w_sum^2;
dN(3,:) = dN(3,:)/w_sum - N*dw_zeta/w_sum^2;

N = N/w_sum;

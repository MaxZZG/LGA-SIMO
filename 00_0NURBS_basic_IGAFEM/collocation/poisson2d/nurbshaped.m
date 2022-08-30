function [N, dN, ddN, coord, J, dxdxi]=nurbshaped(e,u_hat,v_hat,u_knot,v_knot,b_net,p,q,mcp,ncp,element_nod,coord_ij)

%calculate the shape function and second derivatives


nsd = 2; 
% ------------------------------------------------------------------
%     get nurbs coordinates forml local node 1;
%[ien,inn]=genIen_Inn_2D(p,q,mcp,ncp);
ni = coord_ij(element_nod(e,end),1);
nj = coord_ij(element_nod(e,end),2);
%     get u and v coordinates of integration point;
u =((u_knot(ni+1)-u_knot(ni))*u_hat +u_knot(ni+1) + u_knot(ni))/2;
v =((v_knot(nj+1)-v_knot(nj))*v_hat +v_knot(nj+1) + v_knot(nj))/2;

%     evaluate 1d size functions and derivatives each direction;
M=dersbasisfuns(ni,p,mcp,u,u_knot); % calculate in u direction
P=dersbasisfuns(nj,q,ncp,v,v_knot); % calculate in v direction

coord = zeros(nsd,1);
w= zeros((p+1)*(q+1), 1);
cpts = zeros((p+1)*(q+1), 2);
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
     
        w(icount) = b_net(ni+i-p-1,nj+j-q-1,nsd+1);
        cpts(icount, 1) = b_net(ni+i-p-1,nj+j-q-1,1);
        cpts(icount, 2) = b_net(ni+i-p-1,nj+j-q-1,2);
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

N(1,:) = B .* w';
dN(1,:) = dBxi .* w';
dN(2,:) = dBeta .* w';
ddN(1,:) = ddBxi .* w';
ddN(2,:) = ddBxieta .* w';
ddN(3,:) = ddBeta .* w';




% Compute the sums of B-spline functions
w_sum = sum(N(1,:));
dw_xi = sum(dN(1,:));
dw_eta = sum(dN(2,:));
d2w_xi = sum(ddN(1,:));
d2w_xieta = sum(ddN(2,:));
d2w_eta = sum(ddN(3,:));

% Compute NURBS basis functions and its first and second derivatives in
% local coordinates
ddN(1,:) = ddN(1,:)/w_sum - (2*dN(1,:)*dw_xi + N*d2w_xi)/w_sum^2 + 2*N*dw_xi^2/w_sum^3;
ddN(2,:) = ddN(2,:)/w_sum - (dN(1,:)*dw_eta + dN(2,:)*dw_xi + N*d2w_xieta)/w_sum^2 +...
    2*N*dw_xi*dw_eta/w_sum^3;
ddN(3,:) = ddN(3,:)/w_sum - (2*dN(2,:)*dw_eta + N*d2w_eta)/w_sum^2 + 2*N*dw_eta^2/w_sum^3;
dN(1,:) = dN(1,:)/w_sum - N*dw_xi/w_sum^2;
dN(2,:) = dN(2,:)/w_sum - N*dw_eta/w_sum^2;
N = N/w_sum;


% calculate coordinates in physical space
icount = 0;
for j=1:q+1
    for i=1:p+1
        icount = icount+1;
        coord(1) = coord(1) + N(icount)*cpts(icount, 1);
        coord(2) = coord(2) + N(icount)*cpts(icount, 2);
    end
end




% Compute Jacobian matrix
dxdxi = [dN(1,:)*cpts(:,1), dN(2,:)*cpts(:,1); dN(1,:)*cpts(:,2), dN(2,:)*cpts(:,2)];
J = det(dxdxi);


% Set up the Hessian and the matrix of squared first derivatives
d2xdxi2 = [ddN(1,:)*cpts(:,1), ddN(2,:)*cpts(:,1), ddN(3,:)*cpts(:,1);...
    ddN(1,:)*cpts(:,2), ddN(2,:)*cpts(:,2), ddN(3,:)*cpts(:,2)];

dxdxi2 = [dxdxi(1,1)^2, dxdxi(1,1)*dxdxi(1,2), dxdxi(1,2)^2;...
    2*dxdxi(1,1)*dxdxi(2,1), dxdxi(1,1)*dxdxi(2,2)+dxdxi(1,2)*dxdxi(2,1), 2*dxdxi(1,2)*dxdxi(2,2);...
    dxdxi(2,1)^2, dxdxi(2,1)*dxdxi(2,2) dxdxi(2,2)^2];

% Solve for first derivatives in global coordinates 
dN = dxdxi'\dN;

% Solve for second derivatives in global coordinates 
ddN = dxdxi2'\(ddN - d2xdxi2'*dN);


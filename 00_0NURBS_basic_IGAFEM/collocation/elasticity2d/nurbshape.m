% Subroutine eval_SHAPE.f consumes an element number and the coordinates in the
% parent element of an integration point and returns the vector of all local
% basis functions evaluated at the point and the matrix of gradients for all
% nonzero bais functions with respect to parameters u and v and with respect
% to x and y. 
%
%  June 17, 2003
%  J. Austin Cottrell
%  CES Graduate Student
%  Texas Institute for Computational Engineering Science
%  University of Texas at Austin
%
%  Modify to codes Matlab by :
%  Hung Nguyen Xuan
%
%   Faculty of Mathematics & Informatics, University of Natural Sciences
%   Vietnam   National University?HCM

function [shl,coord, shgradl,shgradg,detj]=nurbshape(e,u,v,u_knot,v_knot,b_net,p,q,mcp,ncp,element_nod,coord_ij)

%global nsd nshl p q mcp ncp;

nsd = 2; 
nshl = (p+1)*(q+1);

shl=zeros(nshl,1); 
shgradl=zeros(nshl,nsd); 
denom_sum=0; 
derv_sum_u=0; 
derv_sum_v=0; 

% ------------------------------------------------------------------
%     get nurbs coordinates forml local node 1;
%[ien,inn]=genIen_Inn_2D(p,q,mcp,ncp);
ni = coord_ij(element_nod(e,nshl),1);
nj = coord_ij(element_nod(e,nshl),2);

%     evaluate 1d size functions and derivatives each direction;
order_deriv = 1;
M=dersbasisfuns(ni,p,order_deriv,u,u_knot); % calculate in u direction
N=dersbasisfuns(nj,q,order_deriv,v,v_knot); % calculate in v direction

% form basis functions and derivatives dr./du and dr./dv;
icount = 0;

 for j = 0:q
    for i = 0:p
        icount = icount+1;
%      basis functions;
        shl(icount,1) = M(1,p+1-i)*N(1,q+1-j)*b_net(ni-i,nj-j,nsd+1);
        denom_sum = denom_sum + shl(icount);
%      derivatives;
        shgradl(icount,1) = M(2,p+1-i)*N(1,q+1-j)*b_net(ni-i,nj-j,nsd+1); %u
        derv_sum_u = derv_sum_u + shgradl(icount,1);
        shgradl(icount,2) = M(1,p+1-i)*N(2,q+1-j)*b_net(ni-i,nj-j,nsd+1); %v
        derv_sum_v = derv_sum_v + shgradl(icount,2);
    end
end 

% divide through by denominator;
shgradl(:,1) = shgradl(:,1)/denom_sum -(shl(:)*derv_sum_u)/(denom_sum^2);
shgradl(:,2) = shgradl(:,2)/denom_sum -(shl(:)*derv_sum_v)/(denom_sum^2);
shl = shl/denom_sum;

%  now calculate gradients.;
%    calculate dx/dxi;

dxdxi=zeros(nsd,nsd); 
coord=zeros(nsd, 1);
icount = 0;
 for j = 0: q
    for i = 0: p
        icount = icount + 1;
        coord(1) = coord(1) + b_net(ni-i,nj-j,1)*shl(icount);
        coord(2) = coord(2) + b_net(ni-i,nj-j,2)*shl(icount);
        dxdxi(1,1) = dxdxi(1,1) + b_net(ni-i,nj-j,1)*shgradl(icount,1);
        dxdxi(1,2) = dxdxi(1,2) + b_net(ni-i,nj-j,1)*shgradl(icount,2);
        dxdxi(2,1) = dxdxi(2,1) + b_net(ni-i,nj-j,2)*shgradl(icount,1);
        dxdxi(2,2) = dxdxi(2,2) + b_net(ni-i,nj-j,2)*shgradl(icount,2);
    end
 end

% compute the inverse of deformation gradient and gradient of shapes in physical coordinates;
%dxdxi
%dxidx=pinv(dxdxi);
%shgradg = shgradl*dxidx;
warning('off', 'MATLAB:singularMatrix')
shgradg = shgradl/dxdxi;
%shgradg = shgradl*pinv(dxdxi);
warning('on', 'MATLAB:singularMatrix')
% Note that DetJ resides in common
detj = det(dxdxi);
%  if(detj < 0) % chu y doan code doi dau detj nay khi can thiet.
%  detj = -detj;
%  end
%pause
clear icount
return;
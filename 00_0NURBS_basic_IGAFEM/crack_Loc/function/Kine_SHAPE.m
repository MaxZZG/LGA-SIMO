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
%   Vietnam   National University–HCM

function [shl,shgradl,shgradg,detj]=Kine_Shape(e,u_hat,v_hat,u_knot,v_knot,b_net,p1,q1,ien,inn)
global nsd nshl p q mcp ncp;

shl=zeros(size(ien,2),1); 
shgradl=zeros(size(ien,2),nsd); 
denom_sum=0; 
derv_sum_u=0; 
derv_sum_v=0; 

% ------------------------------------------------------------------
%     get nurbs coordinates forml local node 1;

ni = inn(ien(e,1),1)-(p-p1);
nj = inn(ien(e,1),2)-(q-q1);
%     get u and v coordinates of integration point;
u =((u_knot(ni+1)-u_knot(ni))*u_hat +u_knot(ni+1) + u_knot(ni))/2;
v =((v_knot(nj+1)-v_knot(nj))*v_hat +v_knot(nj+1) + v_knot(nj))/2;

%     evaluate 1d size functions and derivatives each direction;
M=dersbasisfuns(ni,p1,mcp,u,u_knot); % calculate in u direction
N=dersbasisfuns(nj,q1,ncp,v,v_knot); % calculate in v direction

% form basis functions and derivatives dr./du and dr./dv;
icount = 0;

 for j = 0:q1
    for i = 0:p1
        icount = icount+1;
%      basis functions;
        shl(icount,1) = M(1,p1+1-i)*N(1,q1+1-j)*b_net(ni-i,nj-j,nsd+1);
        denom_sum = denom_sum + shl(icount);
%      derivatives;
        shgradl(icount,1) = M(2,p1+1-i)*N(1,q1+1-j)*b_net(ni-i,nj-j,nsd+1); %u
        derv_sum_u = derv_sum_u + shgradl(icount,1);
        shgradl(icount,2) = M(1,p1+1-i)*N(2,q1+1-j)*b_net(ni-i,nj-j,nsd+1); %v
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
icount = 0;
 for j = 0: q1
    for i = 0: p1
        icount = icount + 1;
        dxdxi(1,1) = dxdxi(1,1) + b_net(ni-i,nj-j,1)*shgradl(icount,1);
        dxdxi(1,2) = dxdxi(1,2) + b_net(ni-i,nj-j,1)*shgradl(icount,2);
        dxdxi(2,1) = dxdxi(2,1) + b_net(ni-i,nj-j,2)*shgradl(icount,1);
        dxdxi(2,2) = dxdxi(2,2) + b_net(ni-i,nj-j,2)*shgradl(icount,2);
    end
 end

% compute the inverse of deformation gradient and gradient of shapes in physical coordinates;
dxidx=inv(dxdxi);
shgradg = shgradl*dxidx;
% Note that DetJ resides in common
detj = det(dxdxi);
%  if(detj < 0) % chu y doan code doi dau detj nay khi can thiet.
%  detj = -detj;
%  end
clear icount
return;
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
%   Vietnam   National University HCM

function [R,shgradl,shgradg,shgradg2d,detj]=Kine_Shape_2nd(e,u_hat,v_hat,u_knot,v_knot,b_net,p,q,nshl,mcp,ncp)
global nsd 

shl=zeros(nshl,1);
shgradl=zeros(nshl,nsd);
shgradl2=zeros(nshl,nsd+1);
denom_sum=0;
derv_sum_u=0;
derv_sum_v=0;
derv_sum_uu = 0 ;
derv_sum_vv = 0 ;
derv_sum_uv = 0 ;

% ------------------------------------------------------------------
%     get nurbs coordinates forml local node 1;
[ien,inn]=genIEN_INN_2D_repeated_Ck(p,q,mcp,ncp);
ni = inn(ien(e,1),1);
nj = inn(ien(e,1),2);
%     get u and v coordinates of integration point;
u =((u_knot(ni+1)-u_knot(ni))*u_hat +u_knot(ni+1) + u_knot(ni))/2;
v =((v_knot(nj+1)-v_knot(nj))*v_hat +v_knot(nj+1) + v_knot(nj))/2;

% evaluate 1d size functions and derivatives each direction
% size and M and N = number of derviatives+shape functions, degree of
% poynomial.
% row 1 of M and N => shape functions
% i^{th} row (i > 1) => i^{th} derivative of the shape function
% calculate in u direction
M = dersbasisfuns(ni,p,mcp,u,u_knot) ;

% calculate in v direction
N = dersbasisfuns(nj,q,ncp,v,v_knot) ;

% form basis functions and derivatives dr./du and dr./dv;
icount = 0;

for j = 0:q
    for i = 0:p
        icount = icount+1;
        
        % basis functions
        shl(icount,1) = M(1,p+1-i)*N(1,q+1-j)*b_net(ni-i,nj-j,nsd+1);
        denom_sum = denom_sum + shl(icount);
        
        % first derivatives
        shgradl(icount,1) = M(2,p+1-i)*N(1,q+1-j)*b_net(ni-i,nj-j,nsd+1); %u
        derv_sum_u = derv_sum_u + shgradl(icount,1);
        shgradl(icount,2) = M(1,p+1-i)*N(2,q+1-j)*b_net(ni-i,nj-j,nsd+1); %v
        derv_sum_v = derv_sum_v + shgradl(icount,2);
        
        % second derivatives
        
        % wrt u
        shgradl2(icount,1) = M(3,p+1-i)*N(1,q+1-j)*b_net(ni-i,nj-j,nsd+1) ;
        derv_sum_uu = derv_sum_uu + shgradl2(icount,1) ;
        
        % wrt v
        shgradl2(icount,2) = M(1,p+1-i)*N(3,q+1-j)*b_net(ni-i,nj-j,nsd+1) ;
        derv_sum_vv = derv_sum_vv + shgradl2(icount,2) ;
        
        % cross derivative...wrt uv
        shgradl2(icount,3) = M(2,p+1-i)*N(2,q+1-j)*b_net(ni-i,nj-j,nsd+1) ;
        derv_sum_uv = derv_sum_uv + shgradl2(icount,3) ;
    end
end

% basis functions
R = shl/denom_sum;

% First derivative.....divide through by denominator
shgradl(:,1) = shgradl(:,1)/denom_sum -(shl(:)*derv_sum_u)/(denom_sum^2);
shgradl(:,2) = shgradl(:,2)/denom_sum -(shl(:)*derv_sum_v)/(denom_sum^2);

% % Second derivative....
% % wrt u
% shgradl2(:,1) = (shgradl2(:,1) - 2*derv_sum_u*shgradl(:,1) - ...
%     derv_sum_uu*shl(:))/denom_sum ;
% 
% % wrt v
% shgradl2(:,3) = (shgradl2(:,3) - 2*derv_sum_v*shgradl(:,2) - ...
%     derv_sum_vv*shl(:))/denom_sum ;
% 
% % wrtuv
% shgradl2(:,2) = (shgradl2(:,2) - derv_sum_uv*shl(:) - ...
%     derv_sum_u*shgradl(:,2) - derv_sum_v*shgradl(:,1))/denom_sum ;

% Second derivative....
% wrt u
shgradl2(:,1) = shgradl2(:,1)/denom_sum - 2*derv_sum_u*shgradl(:,1)/denom_sum^2 ...
   - derv_sum_uu*shl(:)/denom_sum^2+ 2*derv_sum_u^2*shl(:)/denom_sum^3 ;

% wrt v
shgradl2(:,2) = shgradl2(:,2)/denom_sum - 2*derv_sum_v*shgradl(:,2)/denom_sum^2 ...
   - derv_sum_vv*shl(:)/denom_sum^2+ 2*derv_sum_v^2*shl(:)/denom_sum^3 ;

% wrtuv
shgradl2(:,3) = shgradl2(:,3)/denom_sum - shgradl(:,1)*derv_sum_v/denom_sum^2 ...
   - shgradl(:,2)*derv_sum_u/denom_sum^2 ... 
   - shl(:)*derv_sum_uv/denom_sum^2+ 2*shl(:)*derv_sum_u*derv_sum_v/denom_sum^3;

% now calculate gradients.;
% calculate dx/dxi;

dxdxi=zeros(nsd,nsd);
icount = 0;
for j = 0: q
    for i = 0: p
        icount = icount + 1;
        dxdxi(1,1) = dxdxi(1,1) + b_net(ni-i,nj-j,1)*shgradl(icount,1);
        dxdxi(1,2) = dxdxi(1,2) + b_net(ni-i,nj-j,1)*shgradl(icount,2);
        dxdxi(2,1) = dxdxi(2,1) + b_net(ni-i,nj-j,2)*shgradl(icount,1);
        dxdxi(2,2) = dxdxi(2,2) + b_net(ni-i,nj-j,2)*shgradl(icount,2);
    end
end

% compute the inverse of deformation gradient and gradient of shapes in physical coordinates;
dxidx=inv(dxdxi);
shgradg = shgradl*dxidx ;

% Note that DetJ resides in common
detj = det(dxdxi);
%  if(detj < 0) % chu y doan code doi dau detj nay khi can thiet.
%  detj = -detj;
%  end
clear icount

% for higher order derivatives
icount = 0 ;
jac1 = zeros(3,2) ;
jac2 = zeros(3,3) ;
for j = 0:q
    for i = 0:p
        icount = icount + 1 ;
        
        % term 1
        jac1(1,1) = jac1(1,1) + b_net(ni-i,nj-j,1)*shgradl2(icount,1) ;
        jac1(1,2) = jac1(1,2) + b_net(ni-i,nj-j,2)*shgradl2(icount,1) ;
        jac1(2,1) = jac1(2,1) + b_net(ni-i,nj-j,1)*shgradl2(icount,2) ;
        jac1(2,2) = jac1(2,2) + b_net(ni-i,nj-j,2)*shgradl2(icount,2) ;
        jac1(3,1) = jac1(3,1) + b_net(ni-i,nj-j,1)*shgradl2(icount,3) ;
        jac1(3,2) = jac1(3,2) + b_net(ni-i,nj-j,2)*shgradl2(icount,3) ;
                
    end
end
jac2 = [dxdxi(1,1)^2       dxdxi(2,1)^2          2*dxdxi(1,1)*dxdxi(2,1);...
        dxdxi(1,2)^2       dxdxi(2,2)^2          2*dxdxi(1,2)*dxdxi(2,2);...
  dxdxi(1,1)*dxdxi(1,2) dxdxi(2,1)*dxdxi(2,2) (dxdxi(2,1)*dxdxi(1,2)+dxdxi(1,1)*dxdxi(2,2))] ;

% jac2 = [dxdxi(1,1)^2 2*dxdxi(1,1)*dxdxi(1,2) dxdxi(1,2)^2;...
%     dxdxi(1,1)*dxdxi(2,1)  (dxdxi(1,1)*dxdxi(2,2)+dxdxi(1,2)*dxdxi(2,1)),...
%     dxdxi(1,2)*dxdxi(2,2);...
%     dxdxi(2,1)^2 2*dxdxi(2,2)*dxdxi(2,1) dxdxi(2,2)^2] ;

  
term1 = shgradg*jac1' ;
term2 = shgradl2 - term1 ;

shgradg2 = inv(jac2)*term2' ;
shgradg2d = shgradg2' ;

end
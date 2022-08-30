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

function [R,dNdxi dNdeta]=NURBS_Shape_2D_Ck(e,u_hat,v_hat,u_knot,v_knot,p,q,b_net)

nshl = (p+1)*(q+1);    % number of local shape functions
shl=zeros(nshl,1);
shgradl=zeros(nshl,2);

denom_sum=0;
derv_sum_u=0;
derv_sum_v=0;


mcp=length(u_knot)-p-1;
ncp=length(v_knot)-q-1;

% ------------------------------------------------------------------
%     get nurbs coordinates forml local node 1;
[ien,inn]=genIEN_INN_2D_repeated_Ck(p,q,mcp,ncp);
ni = inn(ien(e,1),1);
nj = inn(ien(e,1),2);
%     get u and v coordinates of integration point;
u =( (u_knot(ni+1)-u_knot(ni))*u_hat +u_knot(ni+1) + u_knot(ni))/2;
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
        shl(icount,1) = M(1,p+1-i)*N(1,q+1-j)*b_net(ni-i,nj-j);
        denom_sum = denom_sum + shl(icount);
        
        % first derivatives
        shgradl(icount,1) = M(2,p+1-i)*N(1,q+1-j)*b_net(ni-i,nj-j); %u
        derv_sum_u = derv_sum_u + shgradl(icount,1);
        shgradl(icount,2) = M(1,p+1-i)*N(2,q+1-j)*b_net(ni-i,nj-j); %v
        derv_sum_v = derv_sum_v + shgradl(icount,2);       
    end
end

% basis functions
shl = shl/denom_sum;

% First derivative.....divide through by denominator
dRdxi  = shgradl(:,1)/denom_sum -(shl(:)*derv_sum_u)/(denom_sum^2);
dRdeta = shgradl(:,2)/denom_sum -(shl(:)*derv_sum_v)/(denom_sum^2);

for i=1:length(shl)
    R(i)     = shl(length(shl)+1-i);
    dNdxi(i) = dRdxi(length(shl)+1-i);
    dNdeta(i)= dRdeta(length(shl)+1-i);
end

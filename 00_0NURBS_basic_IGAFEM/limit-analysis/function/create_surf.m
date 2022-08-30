function [X,Y,Z] = create_surf(p,q,U,V,CP)
% draws a NURBS surface grid 50x50 lines

mu = length(U);
mv = length(V);
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

check_input(p,mu,nu,q,mv,nv)

eps=10e-10;

l=1;  % counting index of lines
r=(V(mv-q)-V(q+1))/50;    %incremental step for v
v=V(q+1);
while v <= V(mv-q)+eps
  j = findspan1(v,V,nv);
  s=(U(mu-p)-U(p+1))/50;  %incremental step for u
  u=U(p+1);
  k=1;
  while u <= U(mu-p)+eps
    i = findspan1(u,U,nu);
    XYZ(k,l,1:3) = get_point_surf(p,i,u,U,q,j,v,V,CP);
    k=k+1;
    u=u+s;
  end
  l=l+1;
  v=v+r;
end
X = XYZ(:,:,1);
Y = XYZ(:,:,2);
Z = XYZ(:,:,3);


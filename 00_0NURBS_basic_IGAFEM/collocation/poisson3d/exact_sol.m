function [ u, dudx, dudy, dudz, f, alpha, ddudx, ddudy, ddudz, ddudxdy, ddudxdz, ddudydz ] =exact_sol(x,y,z )

%exact solution for the problem -Delta u + alpha*u = f

alpha = 0;
r=1;
R=10;
u = (x.^2+y.^2+z.^2-(R-r)^2).*(x.^2+y.^2+z.^2-(R+r)^2).*sin(x).*sin(y).*sin(z);
dudx = 2*x.*(x.^2+y.^2+z.^2-(R+r)^2).*sin(x).*sin(y).*sin(z)+(2*(x.^2+y.^2+z.^2-(R-r)^2)).*x.*sin(x).*sin(y).*sin(z)+...
    (x.^2+y.^2+z.^2-(R-r)^2).*(x.^2+y.^2+z.^2-(R+r)^2).*cos(x).*sin(y).*sin(z);
dudy = 2*y.*(x.^2+y.^2+z.^2-(R+r)^2).*sin(x).*sin(y).*sin(z)+(2*(x.^2+y.^2+z.^2-(R-r)^2)).*y.*sin(x).*sin(y).*sin(z)+...
    (x.^2+y.^2+z.^2-(R-r)^2).*(x.^2+y.^2+z.^2-(R+r)^2).*sin(x).*cos(y).*sin(z);
dudz = 2*z.*(x.^2+y.^2+z.^2-(R+r)^2).*sin(x).*sin(y).*sin(z)+(2*(x.^2+y.^2+z.^2-(R-r)^2)).*z.*sin(x).*sin(y).*sin(z)+...
    (x.^2+y.^2+z.^2-(R-r)^2).*(x.^2+y.^2+z.^2-(R+r)^2).*sin(x).*sin(y).*cos(z);

ddudx = (2*(x.^2+y.^2+z.^2-(R+r)^2)).*sin(x).*sin(y).*sin(z)+8*x.^2.*sin(x).*sin(y).*sin(z)+4*x.*(x.^2+y.^2+z.^2-(R+r)^2).*cos(x).*sin(y).*sin(z)+...
    (2*(x.^2+y.^2+z.^2-(R-r)^2)).*sin(x).*sin(y).*sin(z)+(4*(x.^2+y.^2+z.^2-(R-r)^2)).*x.*cos(x).*sin(y).*sin(z)-(x.^2+y.^2+z.^2-(R-r)^2).*(x.^2+y.^2+z.^2-(R+r)^2).*sin(x).*sin(y).*sin(z);
ddudy = (2*(x.^2+y.^2+z.^2-(R+r)^2)).*sin(x).*sin(y).*sin(z)+8*sin(y).*sin(z).*sin(x).*y.^2+4*y.*(x.^2+y.^2+z.^2-(R+r)^2).*sin(x).*cos(y).*sin(z)+...
    (2*(x.^2+y.^2+z.^2-(R-r)^2)).*sin(x).*sin(y).*sin(z)+(4*(x.^2+y.^2+z.^2-(R-r)^2)).*y.*sin(x).*cos(y).*sin(z)-(x.^2+y.^2+z.^2-(R-r)^2).*(x.^2+y.^2+z.^2-(R+r)^2).*sin(x).*sin(y).*sin(z);
ddudz = (2*(x.^2+y.^2+z.^2-(R+r)^2)).*sin(x).*sin(y).*sin(z)+8*sin(y).*sin(z).*sin(x).*z.^2+4*z.*(x.^2+y.^2+z.^2-(R+r)^2).*sin(x).*sin(y).*cos(z)+...
    (2*(x.^2+y.^2+z.^2-(R-r)^2)).*sin(x).*sin(y).*sin(z)+(4*(x.^2+y.^2+z.^2-(R-r)^2)).*z.*sin(x).*sin(y).*cos(z)-(x.^2+y.^2+z.^2-(R-r)^2).*(x.^2+y.^2+z.^2-(R+r)^2).*sin(x).*sin(y).*sin(z);

ddudxdy = 8*x.*y.*sin(x).*sin(y).*sin(z)+2*x.*(x.^2+y.^2+z.^2-(R+r)^2).*sin(x).*cos(y).*sin(z)+(2*(x.^2+y.^2+z.^2-(R-r)^2)).*x.*sin(x).*cos(y).*sin(z)+...
    2*y.*(x.^2+y.^2+z.^2-(R+r)^2).*cos(x).*sin(y).*sin(z)+(2*(x.^2+y.^2+z.^2-(R-r)^2)).*y.*cos(x).*sin(y).*sin(z)+(x.^2+y.^2+z.^2-(R-r)^2).*(x.^2+y.^2+z.^2-(R+r)^2).*cos(x).*cos(y).*sin(z);

ddudxdz = 8*x.*z.*sin(x).*sin(y).*sin(z)+2*x.*(x.^2+y.^2+z.^2-(R+r)^2).*sin(x).*sin(y).*cos(z)+(2*(x.^2+y.^2+z.^2-(R-r)^2)).*x.*sin(x).*sin(y).*cos(z)+...
    2*z.*(x.^2+y.^2+z.^2-(R+r)^2).*cos(x).*sin(y).*sin(z)+(2*(x.^2+y.^2+z.^2-(R-r)^2)).*z.*cos(x).*sin(y).*sin(z)+(x.^2+y.^2+z.^2-(R-r)^2).*(x.^2+y.^2+z.^2-(R+r)^2).*cos(x).*sin(y).*cos(z);

ddudydz = 8*y.*z.*sin(x).*sin(y).*sin(z)+2*y.*(x.^2+y.^2+z.^2-(R+r)^2).*sin(x).*sin(y).*cos(z)+(2*(x.^2+y.^2+z.^2-(R-r)^2)).*y.*sin(x).*sin(y).*cos(z)+...
    2*z.*(x.^2+y.^2+z.^2-(R+r)^2).*sin(x).*cos(y).*sin(z)+(2*(x.^2+y.^2+z.^2-(R-r)^2)).*z.*sin(x).*cos(y).*sin(z)+(x.^2+y.^2+z.^2-(R-r)^2).*(x.^2+y.^2+z.^2-(R+r)^2).*sin(x).*cos(y).*cos(z);

f = -(ddudx + ddudy + ddudz) + alpha*u;



%f = (1+12*pi^2)*u;
function [ u, dudx, dudy, f, alpha, ddudx, ddudxdy, ddudy] =exact_sol(x,y )

%exact solution for the problem -Delta u + alpha_u = f

alpha = 1;
u = (x.^2+y.^2-1).*(x.^2+y.^2-16).*sin(x).*sin(y);
dudx = sin(y).*(2.*x.*(2.*x.^2+2*y.^2-17).*sin(x)+(x.^4+x.^2.*(2*y.^2-17)+y.^4-17*y.^2+16).*cos(x));
dudy = sin(x).*((16+x.^4-17*y.^2+y.^4+x.^2.*(-17+2*y.^2)).*cos(y)+2*y.*(-17+2*x.^2+2*y.^2).*sin(y));
ddudx = sin(y).*(-((x.^4+x.^2.*(2*y.^2-29)+y.^4-21*y.^2+50).*sin(x)-4*x.*(2*x.^2+2*y.^2-17).*cos(x)));
ddudy = -sin(x).*((x.^4+x.^2.*(2*y.^2-21)+y.^4-29*y.^2+50).*sin(y)-4*y.*(2*x.^2+2*y.^2-17).*cos(y));
ddudxdy =  2*x.*sin(x).*((2*x.^2+2*y.^2-17).*cos(y)+4*y.*sin(y))+cos(x).*(2*y.*(2*x.^2+2*y.^2-17).*sin(y)+(x.^4+x.^2.*(2*y.^2-17)+y.^4-17*y.^2+16).*cos(y));

f = -(ddudx + ddudy) + alpha*u;
%f2 = (3*x.^4-67*x.^2-67*y.^2+3*y.^4+6*x.^2.*y.^2+116).*sin(x).*sin(y)+(68*x-8*x.^3-8.*x.*y.^2).*cos(x).*sin(y)+(68*y-8.*y.^3-8*y.*x.^2).*cos(y).*sin(x);


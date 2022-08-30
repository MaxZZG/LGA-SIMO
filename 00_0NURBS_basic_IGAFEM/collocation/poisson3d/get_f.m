function [ f] =get_f(x,y,z )

%exact solution for the problem -Delta u + alpha*u = f

% lenx = length(x);
% leny = length(y);
% lenz = length(z);


% 
% u = x.^2;
% dudx = 2.*x;
% dudy = zeros(leny,lenx);
% f = -2*ones(leny,lenx);

% u = x.*(x-1).*y.*(y-1);
% dudx = (2*x-1).*(y.^2-y);
% dudy = (x.^2-x).*(2*y-1);
% f = -2.*(y.^2-y+x.^2-x);

% u = sin(2*pi*x).*cos(2*pi*y);
% dudx = 2*pi*cos(2*pi*x).*cos(2*pi*y);
% dudy = -2*pi*sin(2*pi*x).*sin(2*pi*y);
% f = 8*pi^2*u;

%alpha = 0;
%u = sin(2*pi*x).*sin(2*pi*y);
%dudx = 2*pi*cos(2*pi*x).*sin(2*pi*y);
%dudy = -2*pi*sin(2*pi*x).*cos(2*pi*y);
%f = 8*pi^2*u;


%u = sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);

%f = (1+12*pi^2)*u;
%f = -2*y.*(y-1).*z.*(z-1)-2.*x.*(x-1).*z.*(z-1)-2*x.*(x-1).*y.*(y-1);

%f = -6*ones(size(x));


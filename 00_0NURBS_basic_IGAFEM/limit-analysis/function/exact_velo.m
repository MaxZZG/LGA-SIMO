function res = exact_velo(x,y)
% res = exact_velo(x,y)
% This function gives exact solution of the proposed problem
% and its derivatives in a given point (x,y)
% res = [u      v
%        u_x    v_x
%        u_y    v_y
%

res = [ 2*x.^2.*(1-x).^2.*y.*(y-1).*(2*y-1)         -2*y.^2.*(1-y).^2.*x.*(x-1).*(2*x-1)
        4*x.*y.*(2*y-1).*(y-1).*(2*x-1).*(x-1)      -2*y.^2.*(y-1).^2.*(6*x.^2-6*x+1)
        2*x.^2.*(6*y.^2-6*y+1).*(x-1).^2            -4*x.*y.*(2*y-1).*(y-1).*(2*x-1).*(x-1) ];
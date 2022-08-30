function [gp,gw]=genGP_GW(ngauss)

[W, Q] = quadrature( ngauss, 'GAUSS', 1 );

gp = Q';
gw = W';

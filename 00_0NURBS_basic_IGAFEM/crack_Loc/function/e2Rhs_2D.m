% Subroutine e3Rhs.m consumes the forces in the x, y, z directions, as well
% as the Gaussian quadrature weight and the shape functions, and assembles
% the local load vector. 
% 
% This is the 3D version
%
% J. Austin Cottrell
% Jurijs Bazilevs
% 
% CAM Graduate Students
% Institute for Computational Engineering Science
% The University of Texas at Austin
%
% Modify to codes Matlab by :
% Hung Nguyen Xuan
%
% Faculty of Mathematics & Informatics, University of Natural Sciences
% Vietnam   National University–HCM


function [rhs,area]=e2Rhs_2D(forx,fory,gwt,detj,shg,rhs,area)
global nshl;
for aa = 1:nshl
    rhs(1,aa) = rhs(1,aa) + forx*shg(aa)*detj*gwt;
    rhs(2,aa) = rhs(2,aa) + fory*shg(aa)*detj*gwt;
end
area = area + detj*gwt;
return






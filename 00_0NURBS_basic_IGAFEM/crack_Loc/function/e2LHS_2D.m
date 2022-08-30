% Subroutine e3LHS.f consumes the parameters lambda and mu, as well as the 
% Gaussian quadrature weight and the spacial gradients, and assembles the
% local stiffness matrix. 
% This is the 3D version
%  
%   J. Austin Cottrell
%   Jurijs Bazilevs
% 
%   CAM Graduate Students
%   Institute for Computational Engineering Science
%   The University of Texas at Austin
%
%  Modify to codes Matlab by :
%  Hung Nguyen Xuan
%
%   Faculty of Mathematics & Informatics, University of Natural Sciences
%   Vietnam   National University–HCM


function xkebe=e2LHS_2D(lambda,mu,gwt,detj,shgradg,xkebe)
global nshl;
% loop over elements in each direction;
for aa = 1:nshl
   for bb = 1: nshl
      xkebe(1,aa,bb) = xkebe(1,aa,bb) +((lambda + 2*mu)*shgradg(aa,1)*shgradg(bb,1) + mu*shgradg(aa,2)*shgradg(bb,2))*detj*gwt;
      xkebe(4,aa,bb) = xkebe(4,aa,bb) +(mu*shgradg(aa,1)*shgradg(bb,1) +(lambda + 2*mu)*shgradg(aa,2)*shgradg(bb,2))*detj*gwt;
      
      xkebe(2,aa,bb) = xkebe(2,aa,bb) +(lambda*shgradg(aa,1)*shgradg(bb,2) +mu*shgradg(aa,2)*shgradg(bb,1))*detj*gwt;
      xkebe(4,aa,bb) = xkebe(4,aa,bb) +(mu*shgradg(aa,1)*shgradg(bb,2) +lambda*shgradg(aa,2)*shgradg(bb,1))*detj*gwt;
   end
end
return
 

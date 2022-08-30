function [rhsb,norsurf]=e2bRhs_2D_pres(shb,nor,detjb,gwt,rhsb,pres,norsurf)
global nshl;
for aa = 1: nshl
    rhsb(1,aa) = rhsb(1,aa) + nor(1)*pres*shb(aa)*detjb*gwt;
    rhsb(2,aa) = rhsb(2,aa) + nor(2)*pres*shb(aa)*detjb*gwt;
end
ndn = sqrt(nor(1)*nor(1)+nor(2)*nor(2));
norsurf = norsurf + ndn*detjb*gwt; % egde length
return;


 
      

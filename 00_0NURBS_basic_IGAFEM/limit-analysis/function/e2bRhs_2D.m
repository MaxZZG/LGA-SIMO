function [rhsb]=e2bRhs_2D(forx,fory,shb,rhsb,detjb,gwt)
global nshl
for aa = 1: nshl
rhsb(1,aa) = rhsb(1,aa) + forx*shb(aa)*detjb*gwt;
rhsb(2,aa) = rhsb(2,aa) + fory*shb(aa)*detjb*gwt;
end
return



function [B, dBdxi]=bernsteinBasis(p,q,xi,eta)
%
% Compute Bernstein basis and first derivaties.
% The basis is defined in [-1,1] to mimic FEM.
%
% Adapted from Ole Jgen Fredheim's master thesis
% VP Nguyen, Cardiff University, UK

%% Initialization
ncpt  = (p+1)*(q+1);
B     = zeros(ncpt,1);
dBdxi = zeros(ncpt,2);

%% Calculation
for j=1:q+1
    for i=1:p+1
        B((p+1)*(j-1)+i)=bernstein(p,i,xi)*bernstein(q,j,eta);
        dBdxi((p+1)*(j-1)+i,1)=0.5*p*(bernstein(p-1,i-1,xi)-bernstein(p- 1,i,xi))*bernstein(q,j,eta);
        dBdxi((p+1)*(j-1)+i,2)=bernstein(p,i,xi)*0.5*q*(bernstein(q-1,j- 1,eta)-bernstein(q-1,j,eta));
    end
end


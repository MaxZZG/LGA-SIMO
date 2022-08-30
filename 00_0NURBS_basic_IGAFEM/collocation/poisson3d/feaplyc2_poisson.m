function [kk,ff]=feaplyc2_poisson(kk,ff,bcdof,bcval)

%----------------------------------------------------------
%  Purpose:
%     Apply constraints to matrix equation [kk]{x}={ff}
%     Assumes homoegeneous boundary conditions
%-----------------------------------------------------------
 
 %v=diag(kk);
 %v(bcdof) = 1;
 kk(bcdof,:) = [];
 kk(:,bcdof) = [];

 ff(bcdof)=[];
 
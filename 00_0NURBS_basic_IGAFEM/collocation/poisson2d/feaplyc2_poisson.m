function [kk,ff]=feaplyc2_poisson(kk,ff,bcdof)

%----------------------------------------------------------
%  Purpose:
%     Apply constraints to matrix equation [kk]{x}={ff}
%     Assume homoegeneous boundary conditions

kk(bcdof,:) = [];
kk(:,bcdof) = [];


ff(bcdof)=[];

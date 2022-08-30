function [kk,ff]=feaplycloc2(kk,ff,bcdof,bcsam,bcvall)

%----------------------------------------------------------
%  Purpose:
%     Apply constraints to matrix equation [kk]{x}={ff} for collocation
%     method, imposing direct dirichlet conditions
%
%  Synopsis:
%     [kk,ff]=feaplybc(kk,ff,bcdof,bcsam,bcval)
%
%  Variable Description:
%     kk - system matrix before applying constraints 
%     ff - system vector before applying constraints
%     bcdof - a vector containing constrained d.o.f
%     bcsam - a vector containing row indices
%     bcval - a vector containing contained value 
%
%     For example, there are constraints at d.o.f=2 and 10
%     and their constrained values are 0.0 and 2.5, 
%     respectively.  Then, bcdof(1)=2 and bcdof(2)=10; and
%     bcval(1)=1.0 and bcval(2)=2.5.
%-----------------------------------------------------------
 

n=length(bcdof);

fac = mean(mean(abs(kk)));
for i=1:n    
    c=bcdof(i);  
    r=bcsam(i);
    kk(r,:) = 0;
    kk(:,c) = 0;
    kk(r,c)=1*fac;
    ff(r)=bcvall(i)*fac;
end



function [kk,mm]=feaplycmass(kk,mm,bcdof)

%  Synopsis:
%     [kk,mm]=feaplybc(kk,mm,bcdof,bcval)
%
%  Variable Description:
%     mm - system matrix before applying constraints 
%     bcdof - a vector containging constrained d.o.f
%-----------------------------------------------------------
 
 n=size(bcdof,2);
 sdof=size(mm,1);

 for i=1:n
    c=bcdof(i);
    for j=1:sdof
       kk(c,j)=0;
       kk(j,c)=0;
       mm(c,j)=0;
       mm(j,c)=0;
    end
	 kk(c,c)=1;
     mm(c,c)=1;
 end


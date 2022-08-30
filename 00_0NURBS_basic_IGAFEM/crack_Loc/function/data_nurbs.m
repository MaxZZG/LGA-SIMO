function [srf]=data_nurbs(Nume,Numr,p,q)
global a b 
mcp=Nume+1;   % u ctrl point
ncp=Numr+1;   % v ctrl point

pnts=zeros(3,mcp,ncp);
wghts=zeros(1,mcp,ncp);

up=linspace(a,b,ncp);   %
for ip=1:length(up)
    [pnt,icrv]=arc(p,up(ip),mcp);
    pnts(:,:,ip)=pnt;
    wghts(:,:,ip)=icrv.coefs(4,:,:);
end
u_knot = [zeros(1,p) linspace(0,1,Nume-p+2) ones(1,p)];
v_knot = [zeros(1,q) linspace(0,1,Numr-q+2) ones(1,q)];
srf = nrbmak(pnts,{u_knot v_knot});
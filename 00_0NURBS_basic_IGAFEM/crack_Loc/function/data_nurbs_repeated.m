function [mcp,ncp,isrf]=data_nurbs_repeated(Nume,Numr,p,q)
global a b order
%         nu = Nume+p-(Nume-1)*(p-1);%nu = 7; 
%         nv = Numr*p-(Numr-1)*(q-1);%nv = 7; 
% k=order;        
        nu = Nume+p-1;
        nv = Numr+q-1;
        
        mcp=nu+1;
        ncp=nv+1;
        
        pnts=zeros(3,mcp,ncp);
        wghts=zeros(1,mcp,ncp);
% up=linspace(R,A,ncp);
% if (mcp<p+1)||(ncp<q+1)
%     error('number of ctrl point must greater than degree');
% end
        up=zeros(1,ncp);
        up(1)=a;
        dd=(b-a)/(nv*ncp/2);
        
for i=2:ncp
    up(i)=up(i-1)+(i-1)*dd;
end

for ip=1:length(up)
    [pnt,icrv]=arc(p,up(ip),mcp);
    pnts(:,:,ip)=pnt;
    wghts(:,:,ip)=icrv.coefs(4,:,:);
end

%-----------------basic mesh-----------------------------------------------
u_knot = [zeros(1,p) linspace(0,1,Nume+1) ones(1,p)];
v_knot = [zeros(1,q) linspace(0,1,Numr+1) ones(1,q)];
srf = nrbmak(pnts,{u_knot v_knot});
% srf.coefs(4,:,:)=wghts;
%-----------------insert loop knots----------------------------------------
ikntu = srf.knots{1}(p+2:length(srf.knots{1})-(p+1));
ikntv = srf.knots{2}(q+2:length(srf.knots{2})-(q+1));
tmp = srf;

    for it=1:p-1
        isrf = nrbkntins(tmp,{ikntu ikntv});
        tmp = isrf;
    end
%-----------------update nurbs properties----------------------------------
mcp = mcp+length(ikntu)*(p-1);
ncp = ncp+length(ikntv)*(q-1);
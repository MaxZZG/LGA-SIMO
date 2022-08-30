function [cpnt,icrv]=arc(p,R,mcp)
%     mcp = nu+1;
    if mcp<p+1
        error('number of ctrl point must greater than degree');
    end
    knt2=[0 0 0 1 1 1]; wm=cos(pi/4);
    x  = R*cos(pi/4);
    y  = R*sin(pi/4);
    xm = x+y*tan(pi/4);
    ctrlpt = [ x wm*xm x;    % w*x - coordinate
               y 0     -y;    % w*y - coordinate
               0 0     0;    % w*z - coordinate
               1 wm    1];   % w   - coordinate
    coefs = zeros(4,3);   % nurbs control points of arc degree 2
    xx = vecrotz(pi/4);
    coefs(:,1:3) = xx*ctrlpt;     % rotate angle
%     xx = vecrotz(pi/2);
    crv = nrbmak(coefs,knt2);   % circular arc degree 2
    crvi = nrbdegelev(crv,p-2);   % elevate degree
    u_knot = [zeros(1,p) linspace(0,1,mcp-p+1) ones(1,p)];
    icrv = nrbkntins(crvi,u_knot(p+2:length(u_knot)-p-1));   % knota addition
%     nrbplots(icrv,20,'curve','b-');
%     plot(icrv.coefs(1,:)./icrv.coefs(4,:),icrv.coefs(2,:)./icrv.coefs(4,:),'ko--')
    cpnt=[];
    for i=1:3
        cpnt=[cpnt; icrv.coefs(i,:)./icrv.coefs(4,:)];
    end
end
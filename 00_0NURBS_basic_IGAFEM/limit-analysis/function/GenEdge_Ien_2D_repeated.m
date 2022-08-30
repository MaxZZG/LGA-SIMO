% This subroutine generates the EDGE_IEN matrix, which relates face 
% numbers and local node numbers to the appropriate global node numbers. 
% It will also generate EDGE_OR which will contain the orientations
% of a given EDGE.
% 
% This version will identify all shape functions with support over 
% an element on the boundary as associated with a particular edge, 
% even though some of these functions will not have support on the
% face itself. This will allow for the calculation of gradients of
% solutions up to the boundary should they be neccessary. This 
% approach should be more robust than merely identifying shape 
% functions on the edge itself. For simplicity, the numbering will 
% not be dependent on edge orientation, as not doing so would 
% introduce a host of different complications in the routines that
% will utilize this connectivity array.
% 
% This is for the 2D code.
% 
% Oct 15, 2009
% 
% Modified to codes Matlab by :
% Hung Nguyen Xuan
%
% Faculty of Mathematics & Informatics, University of Natural Sciences
% Vietnam   National University–HCM

function [edge_ien,edge_or]=GenEdge_Ien_2D(ien,closed_u_flag,closed_v_flag,nedge,numx,numy)
global nshl p q mcp ncp;

edge_ien=zeros(nedge,nshl);
edge_or=zeros(nedge,1);

%    initialize matrices and variables;
edge = 0;

% loop through edge assigning numbers and filling out matrices.;
%           edge numbers will be assigned sequentially starting with the;
%           u edge at v = 1, denoted(u,1). Then we count: the(mcp,v), the(u,ncp), and the(1,v,w),;
%           and lastly(u,v,ocp) surface. once the face has been labled,;
%           we find the number of the element the face belongs to and;
%           assign the local node numbers.;
% 
%        edge orientation scheme:;
%            1 - edge(u,1); eta = c1, canh duoi cung cua square domain
%            2 - edge(mcp,v);xi=c2, canh phai cua square domain in paramater space
%            3 - edge(u,ncp);eta=c3, canh tren cung cua square domain
%            4 - edge(1,v);xi=c4, canh trai cua square domain
%            Chú ý rang, cho bai toan nao co dang tron xuay(diem control point dau va cuoi trung nhau)
%            Vi vay, ta nen ky hieu closed_u_flag=1 de chi rang mat nay la
%            line "dóng" (closed line)
%            Nguoc lai, closed_u_flag=0 is called "opened line"

% edge 1;
% e = 0;
if(closed_v_flag~=1)
  for i = 1:numx
           edge = edge + 1;% edge number
           edge_or(edge) = 1;% edge orientation
           e = i;% element number
%            e = numx*numy+1-i;
           edge_ien(edge,:) = ien(e,:);
  end
end
%  edge 2;
if(closed_u_flag~=1)
 for i = 1:numy
     edge = edge+1;
     edge_or(edge) = 2;
     e =i*numx;
%      e = numx*(numy-1)+1-(i-1)*numx;
     edge_ien(edge,:) = ien(e,:);
  end 
end

%  orientation 3 (eta=ncp);
if(closed_v_flag~=1)
  for i = 1:numx
        edge = edge+1;
        edge_or(edge) = 3;
        e = numx*numy+1-i;
%         e=i;
        edge_ien(edge,:) = ien(e,:);
  end 
end

% orientation 4;(xi=1);
if(closed_u_flag~=1)
   for i = 1:numy
       edge = edge+1;
       edge_or(edge) = 4;
       e = numx*(numy-1)+1-(i-1)*numy;
%        e =i*numx;
       edge_ien(edge,:) = ien(e,:);
    end 
end
return



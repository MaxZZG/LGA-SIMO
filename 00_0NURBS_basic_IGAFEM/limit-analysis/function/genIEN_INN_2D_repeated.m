% This subroutine generates the IEN matrix, which relates element numbers
% and local node numbers to the appropriate global node numbers. The
% routine also generates the INN matrix, which relates global node
% number to the "NURBS coordinates" of the node. This routine is for the 3D code
% 
% J. Austin Cottrell
%
% CAM Graduate Student Institute for Computational Engineering Science
% The University of Texas at Austin
% 
% Modify to codes Matlab by :
% Hung Nguyen Xuan
%
% Faculty of Mathematics & Informatics, University of Natural Sciences
% Vietnam   National University–HCM

function [ien,inn]=genIEN_INN_2D_repeated(p,q,mcp,ncp)
%initialize matrices 
g = 0;
e = 0;
%    loop through control points assigning global node;
%    numbers and filling out ien and inn as we go;
 for j = 1:ncp    % loop through control points in V direction
   for i = 1:mcp % loop through control points in U direction
       g = g+1;
       inn(g,1) = i;
       inn(g,2) = j;
       %if(((i >=(p+1))&&(j >=(q+1))))
     if  ((i >=(p+1)) && (j >=(q+1)) && ( (p==1) || ((mod(i,p)==1) && (mod(j,p)==1)) ))  % CHO NAY!!!
           e = e +1;
           for loop1 = 0:q % number of local nodes in each direction equals "order+1" of local B-Spline basic function 
             for loop2 = 0:p
                 gtemp = g-mcp*loop1-loop2;
                 ln =(p+1)*loop1 +loop2 + 1;
                 ien(e,ln) = gtemp;
              end
           end
       end
   end
end




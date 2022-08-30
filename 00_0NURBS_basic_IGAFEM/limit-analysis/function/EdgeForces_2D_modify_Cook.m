% This subroutine calculates the contributions to the left hand side forcing 
% vector due to boundary edges.
%
% This is for the 2D case.     
%      
%
% Modify to codes Matlab by :
% Hung Nguyen Xuan  
%
% Faculty of Mathematics & Informatics, University of Natural Sciences
% Vietnam   National University–HCM


function RHSG=EdgeForces_2D_modify_Cook(ngauss,nedge,inn,edge_ien,edge_or,ibc_edge,ibc,b_net,RHSG)
global nsd nshl u_knot v_knot p;
global E poisson L D P
%will be the shape function array while shbg will hold the gradients of the shape functions

%***************************************************
% loop over faces;

for ifac = 1: nedge
    area_flag = 0; % initialize flag
%     sctr=edge_ien(ifac,:);           % element edge scatter vector
%     nn=length(sctr);
    sctr=edge_ien(ifac,1:p+1:(p+1)*(p+1));          % element edge scatter vector
    nn=(p+1);
    sctrx = 2.*sctr-1 ;
    sctry = 2.*sctr ;
%Check if face is of nonzero area

% get NURBS coordinates
   ni = inn(edge_ien(ifac,1),1);
   nj = inn(edge_ien(ifac,1),2);

   if((edge_or(ifac)==1)||(edge_or(ifac)==3))
    %face_ora=face_or(ifac)
      da =(u_knot(ni+1) - u_knot(ni))/2;% change in gwt due to mapping from [-1,1]

     if((u_knot(ni)==(u_knot(ni+1))))
       area_flag = 1; % indicates face has zero area then
     end
   else%if((edge_or(ifac)==2)||(edge_or(ifac)==4)) 
      da =(v_knot(nj+1) - v_knot(nj))/2; % change in gwt due to mapping from [-1,1]
      if(v_knot(nj)==(v_knot(nj+1)))
         area_flag = 1; % indicates face has zero area then
      end
    end

%  if face has zero area, skip entirely
 if(area_flag==0)
     
%   check form loads applied to the result faces;
   if(((ibc_edge(ifac,1)==2)||(ibc_edge(ifac,2)==2)))
% skip boundary integrals if load not prescr.
% set up element parameters;
% get gauss point/weight arrays;

      [gp,gw]=genGP_GW(ngauss);
      rhsb=zeros(nsd,nshl); 
  
%      loop over integration points;

      for igaussb = 1: ngauss
      % get boundary size functions and set jacobian;
          [shl,shbg,norm_vec,detjb]=Kine_EDGE_2D(ifac,gp(igaussb),ni,nj,u_knot,v_knot,b_net,edge_or);
         
          RA=shl(1,1:p+1:(p+1)*(p+1));
      % get physical coordinates field;
         % xloc=Physical_coord(b_net,ni,nj,shl);
          gwt = gw(igaussb)*da;
      % anatilycal solution for infinitive plate with hole
          %[ue,estress0]=Analytic_sol_platehole(xloc(1),xloc(2),E,poisson,R);
      % Exact solution of cantilever beam
           fyPt=P;          % y traction at quadrature point
          %RHSG(sctry)=RHSG(sctry)+shl'*fyPt*gwt*detjb; %stiffness matrix
          RHSG(sctry)=RHSG(sctry)+RA'*fyPt*gwt*detjb; %stiffness matrix

      end

   end  %check form loads
 end % end if 
end % end for
return










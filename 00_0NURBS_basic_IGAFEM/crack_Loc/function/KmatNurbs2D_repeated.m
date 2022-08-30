function [LHSK]=KmatNurbs2D(ngauss,nel,inn,ien,b_net,C,LHSK)
global nsd nshl u_knot v_knot p q mcp ncp

% get gaussian points and weights;
[gp,gw]=genGP_GW(ngauss);

area = 0;% Volume of the solid (for debugging)
nel_nza = 0;% Elements of non-zero area

%**************************************************************************
tol=1e-8;% in order to check u_knot(ni) matches u_knot(ni+1) or not ...
% loop over elements;
for iel = 1:nel
  
  sctr=ien(iel,:);           % element scatter vector
  nn=length(sctr);
  sctrB(1:2:2*nn) = 2.*sctr-1 ;
  sctrB(2:2:2*nn) = 2.*sctr   ;
  %  check to see ifmlv current element has nonzero area;
  ni = inn(ien(iel,1),1);% get NURBS coordinates
  nj = inn(ien(iel,1),2);
% element has positive area in the parametric domain
  if(abs(u_knot(ni)-u_knot(ni+1))>tol)&&(abs(v_knot(nj)-v_knot(nj+1))>tol)
     nel_nza = nel_nza + 1;
     % used in calculating quadrature points. The factor of 4 comes from mapping from the [-1,1]
     % line onto a real segment...(in jacobian det)
     da =(u_knot(ni+1) - u_knot(ni))*(v_knot(nj+1) - v_knot(nj))/4;
     %  set up element parameters;
     fx = 0; % body force contributions
     fy = 0;
     xkebe=zeros(2*nshl,2*nshl); % initial element stiffness matrix
     
     %--------------------------------------------------------
     % loop over integration points(ngauss in each direction);
     B=zeros(3,2*nn);
     for igauss = 1: ngauss
        for jgauss = 1: ngauss
%             [N,dNdxi,dNdxy,detj]=Kine_SHAPE_2D_repeated(iel,gp(igauss),gp(jgauss),u_knot,v_knot,b_net);
            [N,dNdxi,dNdxy,detj]=Kine_SHAPE_2D_repeated(iel,gp(igauss),gp(jgauss),u_knot,v_knot,b_net,p,q,nshl,mcp,ncp);
            % calculate given element stiffness matrix and force vector;
              gwt = gw(igauss)*gw(jgauss)*da;
            
              B(1,1:2:2*nn-1) = dNdxy(:,1)';
              B(2,2:2:2*nn)   = dNdxy(:,2)';
              B(3,1:2:2*nn-1) = dNdxy(:,2)';
              B(3,2:2:2*nn)  =  dNdxy(:,1)';

              area = area + detj*gwt;
              xkebe=xkebe+B'*C*B*gwt*detj; %stiffness matrix
%              xme=xme+N'*rho*N*gwt*detj;  % mass matrix 
%             xkebe=e2LHS_2D(lambda,mu,gwt,detj,dNdxy,xkebe); % local stiffness
%            [rhs,area]=e2Rhs_2D(fx,fy,gwt,detj,N,rhs,area);
             clear N dNdxi dNdxy detj 
         end
     end % end integration loop
     LHSK(sctrB,sctrB)=LHSK(sctrB,sctrB)+xkebe; %stiffness matrix
     
  end %end if
end % end loop over elements
disp('Area of plane solid:'); AREA=area
return

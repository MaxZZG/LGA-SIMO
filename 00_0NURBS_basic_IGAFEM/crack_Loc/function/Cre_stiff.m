function [Cv]=Cre_convection(Cv,ien, inn, ngauss, U)
global u_knot v_knot b_net p q mcp ncp 
nnelU = (p+1)*(q+1);
nel=size(ien,1);
[gp,gw]=genGP_GW(ngauss);
tol=1e-8;% in order to check u_knot(ni) matches u_knot(ni+1) or not ...
% loop over elements;
for iel = 1:nel  
  sctr=ien(iel,:);           % element scatter vector
  nn=length(sctr);
  edofU(1:2:2*nn) = 2.*sctr-1 ;
  edofU(2:2:2*nn) = 2.*sctr   ;
  %  check to see ifmlv current element has nonzero area;
  ni = inn(ien(iel,1),1);% get NURBS coordinates
  nj = inn(ien(iel,1),2);
% element has positive area in the parametric domain
  if(abs(u_knot(ni)-u_knot(ni+1))>tol)&&(abs(v_knot(nj)-v_knot(nj+1))>tol)
     da =(u_knot(ni+1) - u_knot(ni))*(v_knot(nj+1) - v_knot(nj))/4;
    
     %  set up element parameters;
     for igauss = 1: ngauss
        for jgauss = 1: ngauss
            gwt = gw(igauss)*gw(jgauss)*da; 
            [N,dNdxi,dNdxy,detj]=Kine_SHAPE_2D_repeated(iel,gp(igauss),gp(jgauss),u_knot,v_knot,b_net,p,q,nnelU,mcp,ncp);
                        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            N_gp= [reshape([1;0]*N',1,2*nn); reshape([0;1]*N',1,2*nn)];
            dNx = [reshape([1;0]*dNdxy(:,1)',1,2*nn); reshape([0;1]*dNdxy(:,1)',1,2*nn)];
            dNy = [reshape([1;0]*dNdxy(:,2)',1,2*nn); reshape([0;1]*dNdxy(:,2)',1,2*nn)];
            velo= N'*[U(2.*sctr-1), U(2.*sctr)];     %[u_x, u_y]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
            Cv(edofU,edofU)=Cv(edofU, edofU)+N_gp'*(velo(1)*dNx+velo(2)*dNy)*detj*gwt;
            
         end
     end % end integration loop
  end %end if
end % end loop over elements
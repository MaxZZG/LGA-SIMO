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


function F=EdgeForces_2D(ngauss,nedge,inn,edge_ien,edge_or,ibc_edge,ibc,b_net,F)
global nsd nshl u_knot v_knot;
global E0 nu0 P;

% loop over faces;
norlength = 0;

for ifac = 1: nedge
    area_flag = 0; % initialize flag
    
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
    else %if((edge_or(ifac)==2)||(edge_or(ifac)==4))
        da =(v_knot(nj+1) - v_knot(nj))/2; % change in gwt due to mapping from [-1,1]
        if(v_knot(nj)==(v_knot(nj+1)))
            area_flag = 1; % indicates face has zero area then
        end
    end
    
    % redefine distribute load for this problems
%     if(edge_or(ifac)==1)
%        P=SigmaInf;
%     end
%     if(edge_or(ifac)==3)
%         P=SigmaInf;
%     end
    
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
                % get physical coordinates field;
%                    xloc=Physical_coord(b_net,ni,nj,shl);
                gwt = gw(igaussb)*da;
                % anatilycal solution for infinitive plate with hole
%                 [ue,estress0]=Analyticsol_platehole(xloc(1),xloc(2),E0,nu0,R);
                % Exact solution of cantilever beam
                %    [eldisp,estress0]=analytic_sol_Cantilever(xloc(1),xloc(2),P,E,poisson,D,L);
                %    taux=estress0(1)*norm_vec(1)+ estress0(3)*norm_vec(2);
                %    tauy=estress0(3)*norm_vec(1)+ estress0(2)*norm_vec(2);
                
                % Distribute for membrane
%                 taux=0;
%                 tauy=SigmaInf;
                % assemble local load vector;
%                 rhsb=e2bRhs_2D(taux,tauy,shl,rhsb,detjb,gwt);
                %-----------------------------------------------------------------------c
                % User-Defined Pressurization                                   
                % Set Value of Pres                                        
                Press = -P;
                [rhsb,norlength]=e2bRhs_2D_pres(shl,norm_vec,detjb,gwt,rhsb,Press,norlength);
                clear ue estress0
            end
            
            % check if displacements are prescribed and assemble into globalml force vector;
            for i = 1:nshl
                g = edge_ien(ifac,i); % global node number
                if(ibc(g,1)==1) % in x-direction
                    rhsb(1,i) = 0;
                end
                if(ibc(g,2)==1) % in y-direction
                    rhsb(2,i) = 0;
                end
                % assemble
                F(nsd*g-1,1) = F(nsd*g-1,1) + rhsb(1,i);
                F(nsd*g,1) = F(nsd*g,1) + rhsb(2,i);
            end
        end
    end % end if
end % end for
length_edge = norlength;
return










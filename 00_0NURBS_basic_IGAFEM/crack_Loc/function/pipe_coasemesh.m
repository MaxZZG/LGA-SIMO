function [CP,U,V,p,q]=pipe_coasemesh(R,r,deg,theta)

%%%%%%%% section of annular

% deg=2; ref=2;  %%% refinement
% theta=90*pi/180;
% R=2;
% r=1;
p=2;
q=1;
U=[0 0 0 1 1 1];
V=[0 0 1 1];
%Control Point coordinates
CP(:,:,1)=[R r; R r; R*cos(theta) r*cos(theta)];
CP(:,:,2)=[0 0; R*tan(theta/2) r*tan(theta/2); R*sin(theta) r*sin(theta)];
CP(:,:,3)=[0 R; 0 R ; 0 R];
CP(:,:,4)=[1 1;  cos(theta/2)  cos(theta/2); 1 1];
[CP,U,V,p,q] = degree_elevate_surf_repeated(p,q,U,V,CP,deg-p,deg-q);

%   for ii = 1:size (CP, 1)
%       X= reshape (CP(ii,:,1), 1, []);
%       Y= reshape (CP(ii,:,2), 1, []);
%       Z= reshape (CP(ii,:,3), 1, []);
%       plot3 (X,Y,Z,'k--'); hold on
%       plot3 (X,Y,Z,'ro','MarkerSize',7,'MarkerFaceColor','r')
%   end
%   
%     for jj = 1:size (CP, 2)
%       X= reshape (CP(:,jj,1), 1, []);
%       Y= reshape (CP(:,jj,2), 1, []);
%       Z= reshape (CP(:,jj,3), 1, []);
%       plot3 (X,Y,Z,'k--')
%       plot3 (X,Y,Z,'ro','MarkerSize',7,'MarkerFaceColor','r')
%     end


% R1 = refinement_vec_repeated_p2(U,ref);
% R2 = refinement_vec_repeated_p2(V,ref);
% 
% [CP,u_knot,v_knot] = knot_refine_surf(p,q,U,V,CP,R1,R2);
% 
% plotNURBS_surf_El_CP(p,q,u_knot,v_knot,CP); hold on
% plot_ctrlnet(CP);
% view(2)
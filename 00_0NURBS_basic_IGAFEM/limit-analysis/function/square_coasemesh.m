function [CP,U,V,p,q]=square_coasemesh(a,b,deg)

%%%%%%%% section cylinder
% R=2;
% r=1;
p=1;
q=1;
U=[0 0 1 1];
V=[0 0 1 1];
%Control Point coordinates
CP(:,:,1)=[0 0; a a];
CP(:,:,2)=[0 b; 0 b];
CP(:,:,3)=[1 1 ;1 1];
CP(:,:,4)=[1 1 ;1 1];

%=====================================================================
% REFINE
[CP,U,V,p,q] = degree_elevate_surf_repeated(p,q,U,V,CP,deg-p,deg-q);


% R1 = refinement_vec_repeated_p2(U,ref);
% R2 = refinement_vec_repeated_p2(V,ref);
% 
% [CP,u_knot,v_knot] = knot_refine_surf(p,q,U,V,CP,R1,R2);
% 
% plotNURBS_surf_El_CP(p,q,u_knot,v_knot,CP); hold on
% plot_ctrlnet(CP);
% view(2)
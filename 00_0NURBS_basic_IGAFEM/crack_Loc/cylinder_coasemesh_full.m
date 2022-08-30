function [CP,U,V,p,q]=cylinder_coasemesh_full(R,r,deg)

%%%%%%%% section cylinder
 theta=90*pi/180;
 R=2;
 r=1;
p=2;
q=1;

U=[0 0 0 .25 .25 .5 .5 .75 .75 1 1 1];
%U=[0 0 0 1 1 1];
V=[0 0 1 1];
%Control Point coordinates
% CP(:,:,1)=[R r; R r; R*cos(theta) r*cos(theta)];
% CP(:,:,2)=[0 0; R*tan(theta/2) r*tan(theta/2); R*sin(theta) r*sin(theta)];
% CP(:,:,3)=[0 R; 0 R ; 0 R];
% CP(:,:,4)=[1 1;  cos(theta/2)  cos(theta/2); 1 1];
t=tan(theta/2);

CP(:,:,1)=[R r; R r; 0 0; -R -r;-R -r;-R -r; 0 0; R r; R r]; % x coord
CP(:,:,2)=[0 0; R*t r*t; R r;R*t r*t;0 0;-R*t -r*t;-R -r;-R*t -r*t; 0 0];
CP(:,:,3)=[0 0; 0 0 ; 0 0;0 0; 0 0 ; 0 0;0 0; 0 0 ; 0 0];
CP(:,:,4)=[1 1; cos(theta/2) cos(theta/2); 1 1;cos(theta/2)  cos(theta/2); 1 1;cos(theta/2)  cos(theta/2); 1 1;cos(theta/2)  cos(theta/2); 1 1];
%%
%=====================================================================
% REFINE
% [CP,U,V,p,q] = degree_elevate_surf_repeated(p,q,U,V,CP,deg-p,deg-q);


% R1 = refinement_vec_repeated_p2(U,ref);
% R2 = refinement_vec_repeated_p2(V,ref);
% 
% [CP,u_knot,v_knot] = knot_refine_surf(p,q,U,V,CP,R1,R2);
% 
% plotNURBS_surf_El_CP(p,q,u_knot,v_knot,CP); hold on
% plot_ctrlnet(CP);
% view(2)
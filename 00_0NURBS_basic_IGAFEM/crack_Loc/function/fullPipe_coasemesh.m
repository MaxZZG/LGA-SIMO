% function [CP,U,V,p,q]=pipe_coasemesh(R,r,deg,theta)
clear all
%%%%%%%% section of annular

deg=2; ref=2;  %%% refinement
cont=deg-1;
% theta=90*pi/180;
R=2;
r=1;
p=2;
q=1;
U=[0 0 0 1/4 1/4 2/4 2/4 3/4 3/4 1 1 1];
V=[0 0 1 1];
%Control Point coordinates
CP(:,:,1)=[R r; R r; 0 0; -R -r; -R -r; -R -r;  0 0 ;  R r ; R r];
CP(:,:,2)=[0 0; R r; R r;  R r ;  0  0; -R -r; -R -r; -R -r; 0 0];
CP(:,:,3)=ones(9,2);
CP(:,:,4)=[1 1; 1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2); 1 1;
           1/sqrt(2) 1/sqrt(2); 1 1; 1/sqrt(2) 1/sqrt(2); 1 1];
       

[CP,U,V,p,q] = degree_elevate_surf(p,q,U,V,CP,deg-p,deg-q);
plotNURBS_surf(U,V,CP);hold on
R1 = refinement_vec(U,ref);
R2 = refinement_vec(V,ref);

% 
% [CP,U,V,p,q] = degree_elevate_surf_repeated(p,q,U,V,CP,deg-p,deg-q);
% R1 = refinement_vec_repeated(U,ref,p-cont);
% R2 = refinement_vec_repeated(V,ref,q-cont);

[CP,uKnot,vKnot] = knot_refine_surf(p,q,U,V,CP,R1,R2);
figure
plotNURBS_surf(uKnot,vKnot,CP); hold on
% plotNURBS_surf_El_CP(p,q,uKnot,vKnot,CP); hold on
view(2)


function plotNURBS_surf_El_CP(p,q,U,V,CP)
% plots the surface, elements and control points

[X,Y,Z] = create_surf(p,q,U,V,CP);

% geometry
 surf(X,Y,Z,'FaceColor','none','EdgeColor','none');
 
xlabel('x'); ylabel('y'); zlabel('z');
hold on;

% element edges
create_el_edges(p,q,U,V,CP)

axis equal;

hold off;

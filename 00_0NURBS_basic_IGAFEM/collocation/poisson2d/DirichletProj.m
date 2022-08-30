function [dirichlet_nodes] = DirichletProj(dirichlet, element_nod, numnodes, knotu, knotv, ngaussedge, nument, b_net, p, q, lenu, lenv, coord_ij, coordinates)
%L2 projection on boundary data

size_dirichlet = size(dirichlet, 1);
num_dirichlet_nodes = 0;

%count the number of nodes on the dirichlet boundary (including
%duplicates) by looping throught the elements
for i=1:size_dirichlet
    ort = dirichlet(i,5);    
    if (ort == 1) || (ort == 3)          
        %top and bottom edges have p+1 functions with support on the        
        %dirichlet boundary
        num_dirichlet_nodes = num_dirichlet_nodes + p + 1;        
    else
        % left and right edges 
        num_dirichlet_nodes = num_dirichlet_nodes + q + 1;
    end
end
  
dirichlet_nodes = zeros(1, num_dirichlet_nodes);
node_counter = 0;
for i=1:size_dirichlet
    ort = dirichlet(i,5);
    element_index = dirichlet(i,1);
    if ort == 1
        dirichlet_nodes(node_counter+1:node_counter+p+1) = element_nod(element_index, 1:1:p+1);
        node_counter = node_counter + p + 1;
    elseif ort == 2
        dirichlet_nodes(node_counter+1:node_counter+q+1) = element_nod(element_index, p+1:p+1:(p+1)*(q+1));
        node_counter = node_counter + q + 1;
    elseif ort == 3        
        dirichlet_nodes(node_counter+1:node_counter+p+1) = element_nod(element_index, q*(p+1)+1:1:(q+1)*(p+1));
        node_counter = node_counter + p + 1;
    elseif ort == 4        
        dirichlet_nodes(node_counter+1:node_counter+q+1) = element_nod(element_index, 1:p+1:q*(p+1)+1);
        node_counter = node_counter + q + 1;
    end
end

%remove duplicate entries
dirichlet_nodes = unique(dirichlet_nodes);

% num_dirichlet_nodes = length(dirichlet_nodes);
% 
% %intialize the mass matrix and rhs vector for L2 projection
% 
% mass = zeros(num_dirichlet_nodes, num_dirichlet_nodes);
% rhs = zeros(num_dirichlet_nodes, 1);
% 
% for i = 1:size_dirichlet
%     ort = dirichlet(i,5);
%     element_index = dirichlet(i,1);
%     shapecount = 0;
%     switch ort 
%         case 1
%             ind_nod_list = 1:1:p+1;
%         case 2
%             ind_nod_list = (p+1):(p+1):(p+1)*(q+1);
%         case 3
%             ind_nod_list = q*(p+1)+1:1:(q+1)*(p+1);
%         case 4
%             ind_nod_list = 1:p+1:q*(p+1)+1;
%         otherwise
%             disp('Unknown orientation')
%     end
%         
%     
%     corient = dirichlet(i, 5);
%     curelm = dirichlet(i, 1);
%     curnod = element_nod(curelm, nument);
%     ni = coord_ij(curnod, 1);
%     nj = coord_ij(curnod, 2);
%     
%     
%     if (corient == 1) || (corient == 3)
%         scalefac = (knotu(ni+1) - knotu(ni))/2;
%     else
%         scalefac = (knotv(nj+1) - knotv(nj))/2;
%     end
%     [gp, gw] = genGP_GW(ngaussedge);
%     
%     localmass = zeros(shapecount, shapecount);
%     localrhs = zeros(shapecount, 1);    
%     scrtx = zeros(shapecount, 1);    
% 
%     %loop over Gauss points and compute the integral
%     for igauss = 1:ngaussedge           
%         [R, coord, ~, J] = nurbedge(gp(igauss), ni, nj, knotu, knotv, b_net, p, q, lenu, lenv, corient);        
%         gwt = gw(igauss)*scalefac;        
%         B = zeros(shapecount, 1);        
%         
%         for j=1:shapecount
%             globnum = shapelist(j,1);
%             enrichnum = shapelist(j,3);
%             enrichind = shapelist(j,2);            
%             locnum = shapelist(j,4);
%             cR = R(locnum);
%             
%             %calculate the value of enrichment and derivatives
%             physx = coord(1);
%             physy = coord(2);
%             xcontrolpt = coordinates(globnum,1);
%             ycontrolpt = coordinates(globnum,2);
%             [zen, ~, ~] = enrichment(physx, physy, enrichnum,xcontrolpt,ycontrolpt);
%             cE = zen;
%             
%             na = cR*cE;            
%             scrtx(j) = dirichlet_node_enrglob{globnum}(enrichind);           
%             B(j) = na;                        
%             
%             %evaluate the prescribed displacement            
%             NodeEval = beamu_d([physx, physy], Emod, nu, L, W, P);
%             if dof==1
%                 %x-direction displacement
%                 ex_disp = NodeEval(1);
%             else
%                 %y-direction displacement
%                 ex_disp = NodeEval(2);
%             end
%             %evaluate the value of the shape function at the
%             %gauss point            
%             localrhs(j) = localrhs(j) + ex_disp.*na.*J.*gwt;
%             
%         end       
%         localmass = localmass + B*B'.*J.*gwt;
%         
%     end    
%     
%     mass(scrtx, scrtx) = mass(scrtx, scrtx) + localmass;
%     rhs(scrtx)=rhs(scrtx)+localrhs;    
% end
% 
% %compute the values of dof by solving the linear system 
% %dirvals = pinv(mass)*rhs;
% dirvals = mass\rhs;

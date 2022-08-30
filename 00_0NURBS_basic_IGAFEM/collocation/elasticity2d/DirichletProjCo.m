function [dirichlet_nodes] = DirichletProjCo(dirichlet, element_nod, p, q)
%returns the indices for the Dirichlet nodes

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

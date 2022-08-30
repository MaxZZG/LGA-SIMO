function [dirichlet_nodes, dirvals] = DirichletProj(dirichlet, dof, element_nod, numnodes, knotu, knotv, ngaussedge, nument, b_net, p, q, lenu, lenv, coord_ij, coordinates, E, nu, tx, rad)

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
    if ort == 1  %bottom
        dirichlet_nodes(node_counter+1:node_counter+p+1) = element_nod(element_index, 1:1:p+1);
        node_counter = node_counter + p + 1;
    elseif ort == 2  %right
        dirichlet_nodes(node_counter+1:node_counter+q+1) = element_nod(element_index, p+1:p+1:(p+1)*(q+1));
        node_counter = node_counter + q + 1;
    elseif ort == 3  %top
        dirichlet_nodes(node_counter+1:node_counter+p+1) = element_nod(element_index, q*(p+1)+1:1:(q+1)*(p+1));
        node_counter = node_counter + p + 1;
    elseif ort == 4      %left
        dirichlet_nodes(node_counter+1:node_counter+q+1) = element_nod(element_index, 1:p+1:q*(p+1)+1);
        node_counter = node_counter + q + 1;
    end
end



%remove duplicate entries
dirichlet_nodes = unique(dirichlet_nodes);
num_dirichlet_nodes = length(dirichlet_nodes);
%map global node index to dirichlet node index
dirichlet_node_glob = zeros(numnodes,1);

glob_dir_num = 0;
for ind_node=1:num_dirichlet_nodes
    cur_node = dirichlet_nodes(ind_node);
    glob_dir_num = glob_dir_num+1;
    dirichlet_node_glob(cur_node) = glob_dir_num;        
end

%intialize the mass matrix and rhs vector for L2 projection

mass = zeros(glob_dir_num, glob_dir_num);
rhs = zeros(glob_dir_num, 1);

for i = 1:size_dirichlet
    ort = dirichlet(i,5);
    element_index = dirichlet(i,1);
    
    switch ort
        case 1  %bottom
            ind_nod_list = 1:1:p+1;
        case 2  %right
            ind_nod_list = (p+1):(p+1):(p+1)*(q+1);
        case 3  %top
            ind_nod_list = q*(p+1)+1:1:(q+1)*(p+1);
        case 4  %left
            ind_nod_list = 1:p+1:q*(p+1)+1;
        otherwise
            disp('Unknown orientation')
    end
    
    %reverse ind_nod_list to match format for nurbedge
    rev_ind_nod_list = nument+1-ind_nod_list;
    
    corient = dirichlet(i, 5);
    curelm = dirichlet(i, 1);
    curnod = element_nod(curelm, nument);
    ni = coord_ij(curnod, 1);
    nj = coord_ij(curnod, 2);
    
    
    if (corient == 1) || (corient == 3)
        scalefac = (knotu(ni+1) - knotu(ni))/2;
    else
        scalefac = (knotv(nj+1) - knotv(nj))/2;
    end
    [gp, gw] = genGP_GW(ngaussedge);
    
    shapecount = length(ind_nod_list);
    localmass = zeros(shapecount, shapecount);
    localrhs = zeros(shapecount, 1);
    scrtx = dirichlet_node_glob(element_nod(element_index, ind_nod_list));
    
    
    
    %loop over Gauss points and compute the integral
    for igauss = 1:ngaussedge
        [R, coord, ~, J] = nurbedge(gp(igauss), ni, nj, knotu, knotv, b_net, p, q, corient);
        gwt = gw(igauss)*scalefac;        
        
        physx = coord(1);
        physy = coord(2);
        
        %evaluate the prescribed displacement
        NodeEval = holeu_d([physx, physy], rad, E, nu, tx);
        if dof==1
            %x-direction displacement
            ex_disp = NodeEval(1);
        else
            %y-direction displacement
            ex_disp = NodeEval(2);
        end
        
        B = R(rev_ind_nod_list);
        
        localrhs = localrhs + ex_disp.*B'.*J.*gwt;                        
        localmass = localmass + B'*B.*J.*gwt;
        
    end
    
    mass(scrtx, scrtx) = mass(scrtx, scrtx) + localmass;
    rhs(scrtx)=rhs(scrtx)+localrhs;
end

%compute the values of dof by solving the linear system
dirvals = mass\rhs;
%resid_dir = max(abs(mass*dirvals-rhs))
function elemDisp = element_disp(e,pos,enrich_node,U)

% From the unknowns vector U, extract the parameters
% associated with the element "e"
% Then epsilon = B*U
% Used for enriched IGA
% Vinh Phu Nguyen
% Johns Hopkins University

global controlPts element

sctr = element(e,:);
nn   = length(sctr);

% stdU contains true nodal displacement

idx    = 0 ;
stdU   = zeros(2*nn,1);

for in = 1 : nn
    idx   = idx + 1;
    nodeI = sctr(in) ;
    stdU(2*idx-1) = U(2*nodeI-1);
    stdU(2*idx)   = U(2*nodeI  );
end

% A contains enriched dofs
A = [];

if ( any(enrich_node(sctr)) == 0 ) % Non-enriched elements
    elemDisp = stdU ; return
else                               % having enriched DOFs
    for in = 1 : nn
        nodeI = sctr(in);
        posI  = pos(nodeI);
        
        if     (enrich_node(nodeI) == 1)     % H(x) enriched node
            AA = [U(2*posI-1);U(2*posI)];
            A  = [A;AA];
        elseif (enrich_node(nodeI) == 2)     % B(x) enriched node
            AA = [U(2* posI-1);
                  U(2* posI);
                  U(2*(posI+1)-1);
                  U(2*(posI+1));
                  U(2*(posI+2)-1);
                  U(2*(posI+2));
                  U(2*(posI+3)-1);
                  U(2*(posI+3));];
            A  = [A;AA];
        elseif (enrich_node(nodeI) == 3)     % inclusion enriched node
            AA = [U(2*posI-1);U(2*posI)];
            A  = [A;AA];
        end
    end
end

% total
elemDisp = [stdU;A];

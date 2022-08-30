function sctrB = assembly(e,element,enrich_node,pos)

% determine the scatter vector to assemble K matrix
% e: element under consideration
% element: element connectivity, seems obsolete because it is defined as a
% global variable.
% enrich_node: enrichment type for enriched nodes

sctr = element(e,:);
nn   = length(sctr);

%for k = 1 : nn
%    sctrBfem(2*k-1) = 2*sctr(k)-1 ;
%    sctrBfem(2*k)   = 2*sctr(k)   ;
%end

% equivalent code with above but without loop

sctrBfem           = zeros(1,2*nn);
sctrBfem(1:2:2*nn) = 2*sctr-1;
sctrBfem(2:2:2*nn) = 2*sctr;

% determine the contribution of enriched nodes

enrnodes = enrich_node(sctr);

if ( any(enrnodes) == 0 ) % Non-enriched elements
    sctrB = sctrBfem ;
    clear sctrBfem;
else
    sn = size(find(enrnodes == 1),1); % Heavise enriched
    tn = size(find(enrnodes == 2),1); % tip enriched
    in = size(find(enrnodes == 3),1); % inclusion enriched
    sctrBxfem = zeros(1,2*(sn*1+tn*4+in*1));
    cnt = 0 ;
    
    % loop over nodes of the current element
    
    for k = 1 : nn
        nk   = sctr(k);
        pnk  = pos(nk);
        pnk2 = 2 * pnk;
        
        % encounter a Heaviside enriched or inclusion enriched node
        
        if     ( enrich_node(nk) == 1) || ( enrich_node(nk) == 3)
            cnt = cnt + 1 ;
            sctrBxfem(2*cnt - 1) = pnk2 - 1;
            sctrBxfem(2*cnt    ) = pnk2    ;
            
        % encounter a tip enriched node
            
        elseif ( enrich_node(nk) == 2)
            cnt = cnt + 1 ;
            sctrBxfem(2*cnt - 1) = pnk2 - 1;
            sctrBxfem(2*cnt    ) = pnk2    ;
            
            cnt = cnt + 1 ;
            sctrBxfem(2*cnt - 1) = pnk2 + 1; %2 * (pnk+1) - 1;
            sctrBxfem(2*cnt    ) = pnk2 + 2; %2 * (pnk+1)    ;
            
            cnt = cnt + 1 ;
            sctrBxfem(2*cnt - 1) = pnk2 + 3; %2 * (pnk+2) - 1;
            sctrBxfem(2*cnt    ) = pnk2 + 4; %2 * (pnk+2)    ;
            
            cnt = cnt + 1 ;
            sctrBxfem(2*cnt - 1) = pnk2 + 5; %2 * (pnk+3) - 1;
            sctrBxfem(2*cnt    ) = pnk2 + 6; %2 * (pnk+3)    ;
        end
    end
    sctrB = [ sctrBfem sctrBxfem ];
    
    clear sctrBfem;
    clear sctrBxfem;
end


% If penalty method is use, then one must modify the stiffness matrix and
% the nodal force vector
% $K_{ij}$ = Kij - alpha \int phi_i phi_j d \gamma_u
% fj       = fj  - alpha \int phi_i u_bar d \gamma_u

% Loop over elements along top edge

fu = zeros(noDofs,1);
k  = zeros(noDofs,noDofs);

[ff kk en] = penaltyBoundaryCondition(elRangeU,elConnU,noElemsU,...
    bottomEdgeMesh,bottomPoints,E0,nu0,stressState,...
    sigmato,xTip,seg,cracklength,pos,1);

fu   = fu + ff;
k    = k + kk;

[ff kk en] = penaltyBoundaryCondition(elRangeU,elConnU,noElemsU,...
    topEdgeMesh,topPoints,E0,nu0,stressState,...
    sigmato,xTip,seg,cracklength,pos,1);

fu   = fu + ff;
k    = k + kk;

[ff kk en] = penaltyBoundaryCondition(elRangeV,elConnV,noElemsV,...
    rightEdgeMesh,rightPoints,E0,nu0,stressState,...
    sigmato,xTip,seg,cracklength,pos,1);

fu   = fu + ff;
k    = k + kk;

[W1,Q1] = quadrature(noGPs1, 'GAUSS', 1 );

%% left edge: natural BCs
% Loop over elements along left edge = noElemsV

for e=1:noElemsV
    xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
    conn  = elConnV(e,:);
    sctrx = 2*leftEdgeMesh(e,:)-1;
    sctry = 2*leftEdgeMesh(e,:);
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N dNdxi] = NURBS1DBasisDers(Xi,p,uKnot,weights);
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        x        = N    *leftPoints(conn,:); % global coord of GP
        jacob1   = dNdxi*leftPoints(conn,:);
        J1       = norm (jacob1);
        
        [sigmaxx,sigmayy,sigmaxy] = exactStressModeI(x,E0,nu0,...
            sigmato,xTip,seg,cracklength);
        
        tx = -sigmaxx;
        ty = -sigmaxy;
        
        f(sctrx) = f(sctrx) + N' * tx * J1 * J2 * wt;
        f(sctry) = f(sctry) + N' * ty * J1 * J2 * wt;
    end
end

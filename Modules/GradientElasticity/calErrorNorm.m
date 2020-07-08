function [DNorm, ENorm] = calErrorNorm(Mesh, NURBS, E, U, ExactSolu)
% function [DNorm, ENorm] = calErrorNorm(Mesh, NURBS, E, U, ExactSolu)
% ----------------------------------------------------------------------
% Calculate displacement norm and energy norm for bar problem
% ------------------------------------------------------------------
% Input:
%       Mesh: Mesh structure of the computational domain
%       NURBS: NURBS structure
%       E: Young's modulus
%       A: area of cross section
%       U: nodal displacement vector
% ---------------------------------------------------------------------
% Output:
%       DNorm: displacement error norm
%       ENorm: energy error norm
% --------------------------------------------------------------------

NGPs = NURBS.Order(1) + 1; % number of gauss points

[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order(1), NURBS.NCtrlPts(1), NURBS.KntVect{1}, 1, NGPs, Mesh.NEl);

Weights = NURBS.Weights;
CtrlPts = reshape(NURBS.CtrlPts3D, 3, [])';

dispNorm = 0;
energyNorm = 0;
% Loop over elements (knot spans)
for e = 1 : Mesh.NEl
    % loop over Gauss points
    for qx = 1 : NGPs
        N0 = Nx(e, qx, :, 1); % bspline basis functions
        N1 = Nx(e, qx, :, 2); % 1st derivarive of bspline
        % rationalize to create nurbs basis functions and their derivatives
        [R0, R1] = Rationalize(Weights(Mesh.El(e, :)), N0(:)', N1(:)');
        % gradient of mapping from parameter space to physical space
        dxdxi = R1 * CtrlPts(Mesh.El(e, :), :);
        % compute the jacobian of physical and parameter domain mapping
        J1 = norm(dxdxi);
        % compute derivative of basis functions w.r.t spatial
        % physical coordinates
        dRdx = J1 \ R1;
        % physical coords of gauss points
        Pts = R0 * CtrlPts(Mesh.El(e, :), :);         
        x = Pts(1);   
        % Numerical solution
        strain = dRdx * U(Mesh.El(e, :));                        
        duh = E(x) * strain;
        uh = R0 * U(Mesh.El(e, :));
        % Exact solution
        ue = ExactSolu.displ(x);
        due = ExactSolu.sigma(x);
        
        dispNorm = dispNorm + (ue - uh)' * (ue - uh) * J1 * Jx(e) * Wx(qx);
        energyNorm = energyNorm + (due - duh)' * (due - duh) * J1 * Jx(e) * Wx(qx);
    end
end
DNorm = sqrt(dispNorm);
ENorm = sqrt(energyNorm);
end
function [mcp,ncp,u_knot,v_knot] = repeated_knot(p,q,nx, ny, k);
% p,q : degree of NuRBS basic functions
% nx,ny :: number element in x,y dir
% k:  C^k continuity element
mcp=p+nx;
ncp=q+ny;
% repeated knot value
        u_rep=zeros((p-k)*(mcp-p-1),1)';
        v_rep=zeros((q-k)*(ncp-q-1),1)';
        for i=1:p-k % multiplicity=p
            u_rep(i:p-k:(p-k)*(mcp-p-1))=(1:(mcp-p-1))/(mcp-p);
            v_rep(i:q-k:(q-k)*(ncp-q-1))=(1:(ncp-q-1))/(ncp-q);
        end
        u_knot=[zeros(1,p+1) u_rep ones(1,p+1)]; % vector u_knot
        v_knot=[zeros(1,q+1) v_rep ones(1,q+1)]; % vector v_knot
        
mcp=length(u_knot)-p-1;         % number of control point in u
ncp=length(v_knot)-q-1;         % number of control point in v             

function TsplineExtraction(knot,insert,spans,p,m)

% Max
% Tspline bezier extraction
[knotExtend,nt] = compute_extended_knot_vector(knot,p);

a = p + 1;
b = a + 1;
nb = 1;

C(1,1) = 0;

mbar = p + 2 + nt + m;
ki = 1;
si = 1;


while b < mbar

    


end






end

function [knotExtend,nt] = compute_extended_knot_vector(knot,p)

% Max
front = ones(1,p) * knot(1);
back = ones(1,p) * knot(end);

knotExtend = [front,knot,back];

nt = numel(front) + numel(back);

end













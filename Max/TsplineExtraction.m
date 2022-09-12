
function c = TsplineExtraction(knot,U,spans,p,m)

% Max
% Tspline bezier extraction
[Ubar,nt] = compute_extended_knot_vector(knot,p);

a = p + 1;
b = a + 1;
nb = 1;

c = zeros(4,4);
c(1,nt+1) = 1;

mbar = p + 2 + nt + m;
ki = 1;
si = 1;

while b < mbar
    c(nb+1,:) = 0;
    add = 0;
    
    temp = find(spans==si);
    if si <= m &&  temp(1)== ki % 
        mult = 0;
        add = 1;
    
        Ubar(b+1:mbar+p-m+si) = Ubar(b:mbar+p-m+si-1);
        Ubar(b) = U(si);
        si = si + 1;
    else
        ki = ki + 1;
        i = b;
        while b < m && Ubar(b+1) == Ubar(b)
            b = b + 1;
        end
    end

    if mult < p 
        numer = Ubar(b) - Ubar(a);
        for j = p:-1:mult+1
            alphas(j-mult) = numer / (Ubar(a+j+add)-Ubar(a));
        end
        r = p - mult;
        for j=1:r
            save = r-j+1;
            s = mult + j;
            for k=p+1:-1:s+1
                alpha = alphas(k-s);
                c(nb,k) = alpha*c(nb,k) +  (1-alpha)*c(nb,k-1);
            end
            if b < m
                c(nb+1,save) = c(nb,p+1);
            end
        end
        nb = nb + 1;
        if b < m
            a = b;
            b = b + 1;
        end
    end
end

end

function [knotExtend,nt] = compute_extended_knot_vector(knot,p)

% Max
front = ones(1,p) * knot(1);
back = ones(1,p) * knot(end);

knotExtend = [front,knot,back];

nt = numel(front);

end













function [z1, z2] = f(x, y);
%F   Volume force in considered domain.

lenx = length(x);
leny = length(y);

%zero volume force
z1 = zeros(leny, lenx);
z2 = zeros(leny, lenx);




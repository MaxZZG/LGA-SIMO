function [A2] = added_variables_strain(sdof,k1,k2,B,S)

A2 = sparse(k2,sdof);
count = 1;
for i = 1:S
    Bx = B{i}(1,:); % see how to store this value in evaluated shape function
    By = B{i}(2,:);
    Bxy = B{i}(3,:);
    A2(count,:) = sqrt(3)*(Bx - By)/2;      % ex - ey 
    A2(count+1,:) = sqrt(3)*Bxy/2;        % gxy
    count = count + 2;
end
A2 = sparse([[A2], [sparse(k2,k1)], [-eye(k2)]]);
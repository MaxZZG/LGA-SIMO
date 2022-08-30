function [A2] = added_variables_plainstrain(sdof,k1,k2,B,S)

A2 = sparse(k2,sdof);
count = 1;
for i = 1:S
    Bx = B{i}(1,:); % see how to store this value in evaluated shape function
    By = B{i}(2,:);
    Bxy = B{i}(3,:);
    A2(count,:)   = (Bx + By);      % ex - ey 
    A2(count+1,:) = Bx - By;        % gxy
    A2(count+2,:) = Bxy;          % ex + ey 
    count = count + 3;
end
A2 = sparse([[A2], [sparse(k2,k1)], [-eye(k2)]]);

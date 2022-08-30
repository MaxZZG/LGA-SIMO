function [A2] = added_variables_strain_ex(sdof,k1,k2,B,S)

A2 = sparse(k2,sdof);
count = 1;
for i = 1:S
    Bx = B{i}(1,:); % see how to store this value in evaluated shape function
    By = B{i}(2,:);
    Bxy = B{i}(3,:);
    A2(count,:) = Bx - By;      % ex - ey 
    A2(count+1,:) = Bxy;        % gxy
    A2(count+2,:) = 4*(Bx + By);      % ex + ey 
    
%     A2(count,:) = Bx ;      % ex - ey 
%     A2(count+1,:) = By;        % gxy
%     A2(count+2,:) = 1/sqrt(2)*(Bxy);      % ex + ey 
    count = count + 3;
end
A2 = sparse([[A2], [sparse(k2,k1)], [-eye(k2)]]);

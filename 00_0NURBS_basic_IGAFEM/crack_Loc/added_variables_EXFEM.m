function [A2] = added_variables_EXFEM(sdof,k1,k2,Bstrain,S)

A2 = sparse(k2,sdof);
count = 1;
for i = 1:S
    
    Bx = Bstrain{i}(1,:); % see how to store this value in evaluated shape function
    
    By = Bstrain{i}(2,:);
    Bxy = Bstrain{i}(3,:);
    A2(count,:) = 2/sqrt(3)*Bx + 1/sqrt(3)*By; % 2/sqrt(3)*exx +1/sqrt(3)*eyy
    A2(count+1,:) = By; % eyy only
    A2(count+2,:) = 1/sqrt(3)*Bxy;  % gxy
    count = count + 3;
%     for k = 1:nG % gauss points
%         Ngx = Nex(k,:);
%         Ngy = Ney(k,:);
%         for ind = 1:nd
%             ii = index(ind); % 2*ii-1 = u, 2*ii = v, (u1 v1...un vn)
%             A2(count ,ii*2-1) = 2/sqrt(3)* Ngx(ind); % 2/sqrt(3)*exx +1/sqrt(3)*eyy
%             A2(count ,ii*2) = 1/sqrt(3)* Ngy(ind);
%             A2(count + 1,ii*2) = Ngy(ind); % eyy only
%             A2(count + 2,ii*2-1) = 1/sqrt(3)* Ngy(ind); %gxy
%             A2(count + 2,ii*2) = 1/sqrt(3)* Ngx(ind);
%         end
%         count = count + 3;
%     end
end
A2 = sparse([[A2], [zeros(k2,k1)], [-eye(k2)]]);
%A21 = [[A2], [zeros(k2,k1)], [-eye(k2)]]
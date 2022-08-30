function [bcdof, bcval]=cavityflowbc(lx,ly,ndof) 
bcdof=[];
bcval=[];

for i=2:lx                            %bottom
    bcdof=[bcdof ndof*i-1 ndof*i];
    bcval=[bcval 0 0];
end
for i=1:ly+1                              % left
    bcdof=[bcdof ndof*((i-1)*(lx+1)+1)-1 ndof*((i-1)*(lx+1)+1)];
    bcval=[bcval 0 0];
end

for i=1:ly+1                            % right
    bcdof=[bcdof ndof*((lx+1)*i)-1 ndof*((lx+1)*i)];
    bcval=[bcval 0 0];
end

for i=2:lx                              % top
    bcdof=[bcdof ndof*((lx+1)*ly+i)-1 ndof*((lx+1)*ly+i)];
    bcval=[bcval 1.0 0];
end

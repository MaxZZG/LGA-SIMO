function nodes=connecElementQ4(numx,numy)

%----------------------------------------------------------
%  Purpose:
%     create connectivity of Q4 element
%
%  Synopsis:
%     nodes=connecElementQ4(numx,numy)
%
%  Variable Description:

%     nodes - connective element matrix
%     numx - number of elements in x direction 
%     numy - number of elements in y direction 


for i=1:numy
    for j=1:numx
     nodes(j+numx*(i-1),1)=j+(numx+1)*(i-1);
     nodes(j+numx*(i-1),2)=j+(numx+1)*(i-1)+1;
     nodes(j+numx*(i-1),3)=j+(numx+1)*i+1;
     nodes(j+numx*(i-1),4)=j+(numx+1)*i;

    end
end
function gcoord=meshCookCoord(L1,L2,L3,numx,numy)

for i=1:numx+1
   for j=1:numy+1
        h=L2-(j-1)*(L2-L3)/numy;
        gcoord((numx+1)*(j-1)+i,1)=(i-1)*L1/numx;
        gcoord((numx+1)*(j-1)+i,2)=(i-1)*L1*(h/L1)/numx+(j-1)*L2/numy;
    end
end


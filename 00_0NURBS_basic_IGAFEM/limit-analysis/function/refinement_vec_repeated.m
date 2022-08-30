function Re = refinement_vec_repeated(U,ref,rep)

R=[];
% k = 1;
for i=1:length(U)-1
  if (U(i)~=U(i+1))
    for j = 1:ref-1;
        R(j) = j/ref*(U(i+1)-U(i))+U(i);        
    end 
  end
end
for ii=1:rep
    Re(ii,:)=R;
end
Re=reshape(Re,1,[]);
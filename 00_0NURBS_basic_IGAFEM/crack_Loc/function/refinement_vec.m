function R = refinement_vec(U,ref)

R=[];
k = 1;
for i=1:length(U)-1
  if (U(i)~=U(i+1))
    for j = 1:ref-1;
      R(k) = j/ref*(U(i+1)-U(i))+U(i);
      k = k+1;
    end 
  end
end
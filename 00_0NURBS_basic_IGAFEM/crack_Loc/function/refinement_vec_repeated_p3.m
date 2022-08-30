function R = refinement_vec_repeated_p3(U,ref)

R=[];
k = 1;
for i=1:length(U)-1
  if (U(i)~=U(i+1))
    for j = 1:ref-1;
      R(k) = j/ref*(U(i+1)-U(i))+U(i);
      R(k+1) = j/ref*(U(i+1)-U(i))+U(i);
      R(k+2) = j/ref*(U(i+1)-U(i))+U(i);
      k =k+3;
    end 
  end
end
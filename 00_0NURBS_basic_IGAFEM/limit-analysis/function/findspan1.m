function i = findspan(u,U,n)

m=length(U);
eps=10e-10;
if (abs(u-U(n+1))<eps)  % special case: last knot (open knot vector assumed)
  i = n;
  return
end
for i = 1:(m-1)
  if (u<U(i+1))
    return
  end
end  

error('u outside of U!')
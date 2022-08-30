function check_input(p,mu,nu,q,mv,nv)
% checks compatibility of input parameters

if (nu+p+1 ~= mu)
  error('U, p and Control points dont match!')
end
if (nv+q+1 ~= mv)
  error('V, q and Control points dont match!')
end
function [CPr,Ur,Vr] = knot_refine_surf_repeated(p,q,U,V,CP,Ru,Rv)

% adapted from algorithm in Piegl, Les. "The NURBS Book". Springer-Verlag: 
% Berlin 1995; p. 167.
% Ru, Rv: Vectors with knots to be inserted in U and V

% nu = length(CP(:,1,1));
% nv = length(CP(1,:,1));

nu = length(U)-p-1;
nv = length(V)-q-1;

% form projective control points Pw
for j = 1:nv
  for i = 1:nu
    Pw(i,j,1:3) = CP(i,j,1:3)*CP(i,j,4);
    Pw(i,j,4)   = CP(i,j,4);
  end
end

if isempty(Ru)
  Ur=U;
  Qw=Pw;
else
  r = length(Ru); 
  mu = nu+p+1;
  a = findspan1(Ru(1),U,nu);
  b = findspan1(Ru(r),U,nu)+1;

  for col = 1:nv
    for  j = 1:a-p;    Qw(j,col,:)   = Pw(j,col,:);   end
    for  j = b-1:nu;   Qw(j+r,col,:) = Pw(j,col,:);   end
  end
  for j = 1:a;      Ur(j)   = U(j);  end
  for j = b+p:mu;   Ur(j+r) = U(j);  end
    
  i = b+p-1;   k = i+r;
  for  j = r:-1:1
    while (Ru(j)<=U(i) && i>a)
      for col = 1:nv
        Qw(k-p-1,col,:) = Pw(i-p-1,col,:);
      end
      Ur(k) = U(i);
      k = k-1;   i = i-1;
    end
    
    for col = 1:nv
      Qw(k-p-1,col,:) = Qw(k-p,col,:);
    end
    for l = 1:p
      ind = k-p+l;
      alfa = (Ru(j)-Ur(k+l)) / (U(i-p+l)-Ur(k+l));
      for col = 1:nv
        Qw(ind-1,col,:) = alfa*Qw(ind-1,col,:) + (1-alfa)*Qw(ind,col,:);
      end
    end
    Ur(k) = Ru(j);
    k = k-1;
  end
  Pw = Qw;
  nu = nu+r;
end

if isempty(Rv)
  Vr=V;
else
  r=length(Rv); 
  mv = nv+q+1;
  a = findspan1(Rv(1),V,nv);
  b = findspan1(Rv(r),V,nv)+1;
  
  for row = 1:nu
    for j = 1:a-q;    Qw(row,j,:) = Pw(row,j,:);    end
    for j = b-1:nv;   Qw(row,j+r,:) = Pw(row,j,:);  end
  end 
  for j = 1:a;      Vr(j)   = V(j);  end
  for j = b+q:mv;   Vr(j+r) = V(j);  end

  i = b+q-1;   k = i+r;
  for j = r:-1:1
    while ((Rv(j)<=V(i)) && (i>a))
      for row = 1:nu
        Qw(row,k-q-1,:) = Pw(row,i-q-1,:);
      end
      Vr(k) = V(i);
      k = k-1;   i = i-1;
    end

    for row = 1:nu
      Qw(row,k-q-1,:) = Qw(row,k-q,:);
    end
    for l = 1:q
      ind = k-q+l;
      alfa = (Rv(j)-Vr(k+l)) / (V(i-q+l)-Vr(k+l));
      for row = 1:nu
        Qw(row,ind-1,:) = alfa*Qw(row,ind-1,:) + (1-alfa)*Qw(row,ind,:);
      end
    end
    Vr(k) = Rv(j);
    k = k-1;
  end
end

for j = 1:length(Qw(1,:,1))
  for i = 1:length(Qw(:,1,1))
    CPr(i,j,1:3) = Qw(i,j,1:3)/Qw(i,j,4);
    CPr(i,j,4)   = Qw(i,j,4);
  end
end
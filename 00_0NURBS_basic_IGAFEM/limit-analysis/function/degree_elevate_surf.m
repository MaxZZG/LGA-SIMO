function [CPr,Ur,Vr,pr,qr] = degree_elevate_surf(p,q,U,V,CP,tp,tq)
% elevates polynomial degrees p and q by tp and tq

% adapted from algorithm in Piegl, Les. "The NURBS Book". Springer-Verlag:
% Berlin 1995; p. 206.

nu = length(CP(:,1,1));
nv = length(CP(1,:,1));
for i = 1:nu
  for j = 1:nv
    Pw(i,j,1:3) = CP(i,j,1:3)*CP(i,j,4);
    Pw(i,j,4)   = CP(i,j,4);
  end
end
if (tp>0)
  mu = nu+p+1;
  pr = p+tp;
  ku = 1;
  for i=1:mu-1;  if(U(i)~=U(i+1));  ku = ku+1;   end;   end 
  Ur = zeros(1,mu+tp*ku);

   bezalfs(1,1) = 1.0;
  bezalfs(pr+1,p+1) = 1.0;

  for i = 1:fix(pr/2)   
    inv = 1/nchoosek(pr,i);
    mpi = min(p,i);
    for j = max(0,i-tp):mpi   
      bezalfs(i+1,j+1) = inv*nchoosek(p,j)*nchoosek(tp,i-j);
      bezalfs(pr-i+1,p-j+1) = bezalfs(i+1,j+1);    
    end
  end

  mh = pr;
  kind = pr+1;
  r=-1;
  a=p+1;
  b=p+2;
  cind=1;
  ua = U(a);
  Qw(1,:,:) = Pw(1,:,:);

  for i = 1:pr+1
    Ur(i) = ua;
  end

  for i = 1:p+1
    bpts(i,:,:) = Pw(i,:,:);
  end

  while (b<mu)
    i = b;
    while (b<mu && U(b)==U(b+1))
      b = b+1;      
    end
    mult = b-i+1;   
    mh = mh+mult+tp;
    ub = U(b);
    oldr = r;       
    r = p-mult;     
  
    if (oldr>0);    lbz = fix((oldr+2)/2);
    else           lbz = 1;
    end
  
    if (r>0);       rbz = pr-fix((r+1)/2);
    else           rbz = pr;
    end
  
    alfs = zeros(1,p-mult);
    if (r>0)
      numer = ub - ua;
      for k = p:-1:mult+1
       alfs(k-mult) = numer/(U(a+k)-ua);   
      end
      for j = 1:r   
        save = r-j+1;
        s = mult+j;
        for k = p+1:-1:s+1    
          bpts(k,:,:) = alfs(k-s)*bpts(k,:,:)+(1-alfs(k-s))*bpts(k-1,:,:);
        end
        Nextbpts(save,:,:) = bpts(p+1,:,:);
      end
    end
  
    for i = lbz:pr
      ebpts(i+1,1:nv,1:4) = 0;
      mpi = min(p,i);
      for j = max(0,i-tp):mpi   
        ebpts(i+1,:,:) = ebpts(i+1,:,:) + bezalfs(i+1,j+1)*bpts(j+1,:,:);
      end
    end
  
    if (oldr>1)
      first = kind-2;
      last = kind;
      den = ub-ua;
      bet = (ub-Ur(kind))/den;
      for tr = 1:(oldr-1)
        i = first;
        j = last;
        kj = j-kind+1;
        while ((j-i)>tr)

          if (i<cind)
            alf = (ub-Ur(i+1))/(ua-Ur(i+1));
            Qw(i+1,:,:) = (alf*Qw(i+1,:,:)+(1.0-alf)*Qw(i,:,:));
          end
          if (j>=lbz)
            if ((j-tr)<=(kind-pr+oldr))
              gam = (ub-Ur(j-tr+1))/den;
              ebpts(kj+1,:,:) = gam*ebpts(kj+1,:,:) + (1.0-gam)*ebpts(kj+2,:,:);
            else
              ebpts(kj+1,:,:) = bet*ebpts(kj+1,:,:) + (1.0-bet)*ebpts(kj+2,:,:);
            end
          end
          i = i+1;
          j = j-1;
          kj = kj-1;
        end
        first = first-1;
        last = last+1;
      end
    end   
     
    if (a~=p+1)   
      for i = 0:(pr-oldr-1)
         Ur(kind+1) = ua;
         kind = kind+1;
      end
    end
  
    for j = lbz:rbz   
      Qw(cind+1,:,:) =  ebpts(j+1,:,:);
      cind = cind +1;
    end
  
    if (b<mu)   
      for j = 0:r-1
        bpts(j+1,:,:) = Nextbpts(j+1,:,:);
      end
      for j = r:p
        bpts(j+1,:,:) = Pw(b-p+j,:,:);
      end
      a = b;
      b = b+1;
      ua = ub;
    else   
      for i = 0:pr
        Ur(kind+i+1) = ub;
      end
    end
  end
  
elseif (tp==0)
  Ur=U;  pr=p;
  Qw=Pw;
end

if (tq>0)
  clear ('Pw');
  Pw(:,:,:)=Qw(:,:,:);
  nu = length(Pw(:,1,1));
  clear ('Qw');  
  clear ('bpts');
  clear ('Nextbpts');
  clear ('ebpts');
  clear ('bezalfs');
  
  mv = nv+q+1;
  qr = q+tq;  
  kv = 1;
  for i=1:mv-1;  if(V(i)~=V(i+1));  kv = kv+1;   end;   end 
  Vr = zeros(1,mv+tq*kv);

  bezalfs(1,1) = 1.0;
  bezalfs(qr+1,q+1) = 1.0;

  for i = 1:fix(qr/2)  
    inv = 1/nchoosek(qr,i);
    mpi = min(q,i);
    for j = max(0,i-tq):mpi   
      bezalfs(i+1,j+1) = inv*nchoosek(q,j)*nchoosek(tq,i-j);
      bezalfs(qr-i+1,q-j+1) = bezalfs(i+1,j+1);    
    end
  end

  mh = qr;
  kind = qr+1;
  r=-1;
  a=q+1;
  b=q+2;
  cind=1;
  ua = V(a);
  Qw(:,1,:) = Pw(:,1,:);

  for i = 1:qr+1
    Vr(i) = ua;
  end

  for i = 1:q+1
    bpts(:,i,:) = Pw(:,i,:);
  end

  while (b<mv)
    i = b;
    while (b<mv && V(b)==V(b+1))
      b = b+1;      
    end
    mult = b-i+1;   
    mh = mh+mult+tq;
    ub = V(b);
    oldr = r;       
    r = q-mult;     
  
    if (oldr>0);    lbz = fix((oldr+2)/2);
    else           lbz = 1;
    end
  
    if (r>0);       rbz = qr-fix((r+1)/2);
    else           rbz = qr;
    end
  
    if (r>0)
      numer = ub - ua;
      for k = q:-1:mult+1
       alfs(k-mult) = numer/(V(a+k)-ua);   
      end
      for j = 1:r   
        save = r-j+1;
        s = mult+j;
        for k = q+1:-1:s+1    
          bpts(:,k,:) = alfs(k-s)*bpts(:,k,:)+(1-alfs(k-s))*bpts(:,k-1,:);
        end
        Nextbpts(:,save,:) = bpts(:,q+1,:);
      end
    end
    
    for i = lbz:qr
      ebpts(1:nu,i+1,1:4) = 0;
      mpi = min(q,i);
      for j = max(0,i-tq):mpi   
        ebpts(:,i+1,:) = ebpts(:,i+1,:) + bezalfs(i+1,j+1)*bpts(:,j+1,:);
      end
    end
    
    if (oldr>1)
      first = kind-2;
      last = kind;
      den = ub-ua;
      bet = (ub-Vr(kind))/den;
      for tr = 1:(oldr-1)
        i = first;
        j = last;
        kj = j-kind+1;
        while ((j-i)>tr)
          if (i<cind)
            alf = (ub-Vr(i+1))/(ua-Vr(i+1));
            Qw(:,i+1,:) = (alf*Qw(:,i+1,:)+(1.0-alf)*Qw(:,i,:));
          end
          if (j>=lbz)
            if ((j-tr)<=(kind-qr+oldr))
              gam = (ub-Vr(j-tr+1))/den;
              ebpts(:,kj+1,:) = gam*ebpts(:,kj+1,:) + (1.0-gam)*ebpts(:,kj+2,:);
            else
              ebpts(:,kj+1,:) = bet*ebpts(:,kj+1,:) + (1.0-bet)*ebpts(:,kj+2,:);
            end
          end
          i = i+1;
          j = j-1;
          kj = kj-1;
        end
        first = first-1;
        last = last+1;
      end
    end   
     
    if (a~=q+1)   
      for i = 0:(qr-oldr-1)
         Vr(kind+1) = ua;
         kind = kind+1;
      end
    end
  
    for j = lbz:rbz   
      Qw(:,cind+1,:) =  ebpts(:,j+1,:);
      cind = cind +1;
    end
  
    if (b<mv)  
      for j = 0:r-1
        bpts(:,j+1,:) = Nextbpts(:,j+1,:);
      end
      for j = r:q
        bpts(:,j+1,:) = Pw(:,b-q+j,:);
      end
      a = b;
      b = b+1;
      ua = ub;
    else   
      for i = 0:qr
        Vr(kind+i+1) = ub;
      end
    end
  end

elseif (tq==0)
  Vr=V;  qr=q;
end

for i = 1:length(Qw(:,1,1))
  for j = 1:length(Qw(1,:,1))
    CPr(i,j,1:3) = Qw(i,j,1:3)/Qw(i,j,4);
    CPr(i,j,4)   = Qw(i,j,4);
  end
end

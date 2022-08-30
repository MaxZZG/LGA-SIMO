%script for imposing Dirichlet boundary condtions
%uses interpolation at Control Points


tol = 1e-9;  %tolerance for if statements


%interpolation on the left edge
DirichletNodes = 1:lenu:lenu*(lenv-1)+1;
lendn = length(DirichletNodes);
%globaldir = zeros(dim*lendn,1);

%create an array of Dirichlet nodes in parameter space

%lendn = 2*lendn;
mult_fac = 3; %multiply number of sample points on Dirichlet boundary to improve stability 
DirichletParmNodes = linspace(0,1,mult_fac*lendn);

DirichletOrient = 4; %left edge

%set-up local-global connectivity matrix for interpolation

if DirichletOrient == 4 
    loc2globint = zeros(nument, length(knotv)-1);
    for i=1:length(knotv)-1
        loc2globint(nument:-(p+1):p+1,i)=(i-p):i;
        
    end
end

InterpolationMatrix = zeros(mult_fac*lendn,lendn);
InterpolationRHSX = zeros(mult_fac*lendn, 1);
InterpolationRHSY = zeros(mult_fac*lendn, 1);

for i=1:mult_fac*lendn
    %i
    curnode = DirichletParmNodes(i);
    if DirichletOrient == 4
        knot = knotv;
        ni = p+1;
        for j=1:length(knot)-1
            
            if (abs(knot(j)-knot(j+1))>tol) && ((curnode-knot(j))>-tol) && ((knot(j+1)-curnode) > -tol)
                nj = j;
                %convert to master coordinates on interval [-1, 1]
                coord = (2*curnode-knot(j)-knot(j+1))/(knot(j+1)-knot(j));
                [shb, phyc, normal, detjb] = nurbedge(coord, ni, nj, knotu, knotv, b_net, p, q, DirichletOrient);
                
                for k = p+1:p+1:nument
                    
                    globnum = loc2globint(k, nj);
                    InterpolationMatrix(i,globnum) = shb(k);
                end
                                                
                NodeEval = holeu_d(phyc', rad, Emod, nu, tx);

                InterpolationRHSX(i) = NodeEval(1);
                InterpolationRHSY(i) = NodeEval(2);
            end
        end
    end
end
InterpolationValuesX = InterpolationMatrix\InterpolationRHSX;
InterpolationValuesY = InterpolationMatrix\InterpolationRHSY;

lendn=size(InterpolationValuesX,1);
InterpolationValuesXY = reshape([InterpolationValuesX, InterpolationValuesY]',2*lendn,1);

%--------------------------------------------
%interpolation on the right edge

DirichletNodes = lenu:lenu:lenu*lenv;
lendn = length(DirichletNodes);
%globaldir = zeros(dim*lendn,1);

%create an array of Dirichlet nodes in parameter space

DirichletParmNodes = linspace(0,1,mult_fac*lendn);

DirichletOrient = 2; %right edge

%set-up local-global connectivity matrix for interpolation

if DirichletOrient == 2
    loc2globint = zeros(nument, length(knotv)-1);
    for i=1:length(knotv)-1
      
        loc2globint(nument-p:-(p+1):1,i)=(i-p):i;
    end
end
        
InterpolationMatrix = zeros(mult_fac*lendn,lendn);
InterpolationRHSX = zeros(mult_fac*lendn, 1);
InterpolationRHSY = zeros(mult_fac*lendn, 1);
for i=1:mult_fac*lendn
    %i
    curnode = DirichletParmNodes(i);
    if DirichletOrient == 2
        knot = knotv;
        ni = length(knotu)-p-1;
        for j=1:length(knot)-1

            if abs(knot(j)-knot(j+1))>tol && (knot(j)<=curnode) && (curnode <= knot(j+1))
                nj = j;
                %convert to master coordinates on interval [-1, 1]
                coord = (2*curnode-knot(j)-knot(j+1))/(knot(j+1)-knot(j));
                [shb, phyc, normal, detjb] = nurbedge(coord, ni, nj, knotu, knotv, b_net, p, q, DirichletOrient);
                            
                for k = 1:p+1:nument-p
       
                    globnum = loc2globint(k, nj);
                    InterpolationMatrix(i,globnum) = shb(k);
                end

                NodeEval = holeu_d(phyc', rad, Emod, nu, tx);
                InterpolationRHSX(i) = NodeEval(1);
                InterpolationRHSY(i) = NodeEval(2);
            end
        end
    end
end

                
InterpolationValuesX = InterpolationMatrix\InterpolationRHSX;
InterpolationValuesY = InterpolationMatrix\InterpolationRHSY;

lendn=size(InterpolationValuesX,1);
InterpolationValuesXY = [InterpolationValuesXY; reshape([InterpolationValuesX, InterpolationValuesY]',2*lendn,1)];



%interpolation on the bottom edge
DirichletNodes = 1:1:lenu;
lendn = length(DirichletNodes);


%create an array of Dirichlet nodes in parameter space

DirichletParmNodes = linspace(0,1,mult_fac*lendn);

DirichletOrient = 1; %bottom edge

%set-up local-global connectivity matrix for interpolation

if DirichletOrient == 1
    loc2globint = zeros(nument, length(knotu)-1);
    for i=1:length(knotu)-1
        loc2globint(nument:-1:nument-p,i)=(i-p):i;
    end
end
        
InterpolationMatrix = zeros(mult_fac*lendn,lendn);
InterpolationRHSX = zeros(mult_fac*lendn, 1);
InterpolationRHSY = zeros(mult_fac*lendn, 1);
for i=1:mult_fac*lendn
    %i
    curnode = DirichletParmNodes(i);
    if DirichletOrient == 1
        knot = knotu;
        nj = p+1;
        for j=1:length(knot)-1
            if abs(knot(j)-knot(j+1))>tol && (knot(j)<=curnode) && (curnode <= knot(j+1))
                ni = j;
                %convert to master coordinates on interval [-1, 1]
                coord = (2*curnode-knot(j)-knot(j+1))/(knot(j+1)-knot(j));
                [shb, phyc, normal, detjb] = nurbedge(coord, ni, nj, knotu, knotv, b_net, p, q, DirichletOrient);
                            
                for k = nument-p:1:nument
                    globnum = loc2globint(k, ni);
                    InterpolationMatrix(i,globnum) = shb(k);
                end
                
                NodeEval = holeu_d(phyc', rad, Emod, nu, tx);
                InterpolationRHSX(i) = NodeEval(1);
                InterpolationRHSY(i) = NodeEval(2);
            end
        end
    end
end

                
InterpolationValuesX = InterpolationMatrix\InterpolationRHSX;
InterpolationValuesY = InterpolationMatrix\InterpolationRHSY;

lendn=size(InterpolationValuesX,1);
InterpolationValuesXY = [InterpolationValuesXY; reshape([InterpolationValuesX, InterpolationValuesY]',2*lendn,1)];

%interpolation on the top edge
DirichletNodes = lenu*(lenv-1)+1:1:lenu*lenv;
lendn = length(DirichletNodes);
globaldir = zeros(dim*lendn,1);

%create an array of Dirichlet nodes in parameter space

%lendn = 2*lendn;
DirichletParmNodes = linspace(0,1,mult_fac*lendn);

DirichletOrient = 3; %top edge

%set-up local-global connectivity matrix for interpolation

if DirichletOrient == 3
    loc2globint = zeros(nument, length(knotu)-1);
    for i=1:length(knotu)-1
        loc2globint(p+1:-1:1,i)=(i-p):i;
    end
end
        
InterpolationMatrix = zeros(mult_fac*lendn,lendn);
InterpolationRHSX = zeros(mult_fac*lendn, 1);
InterpolationRHSY = zeros(mult_fac*lendn, 1);
for i=1:mult_fac*lendn
    %i
    curnode = DirichletParmNodes(i);
    if DirichletOrient == 3
        knot = knotu;
        nj = length(knotv)-q-1;
        for j=1:length(knot)-1
            if abs(knot(j)-knot(j+1))>tol && (knot(j)<=curnode) && (curnode <= knot(j+1))
                ni = j;
                %convert to master coordinates on interval [-1, 1]
                coord = (2*curnode-knot(j)-knot(j+1))/(knot(j+1)-knot(j));
                [shb, phyc, normal, detjb] = nurbedge(coord, ni, nj, knotu, knotv, b_net, p, q, DirichletOrient);
               
                for k = 1:1:p+1
                    globnum = loc2globint(k, ni);
                    InterpolationMatrix(i,globnum) = shb(k);
                end
               
                NodeEval = holeu_d(phyc', rad, Emod, nu, tx);
                
                InterpolationRHSX(i) = NodeEval(1);
                InterpolationRHSY(i) = NodeEval(2);
            end
        end
    end
end

                
InterpolationValuesX = InterpolationMatrix\InterpolationRHSX;
InterpolationValuesY = InterpolationMatrix\InterpolationRHSY;

lend=size(InterpolationValuesX,1);
InterpolationValuesXY = [InterpolationValuesXY; reshape([InterpolationValuesX, InterpolationValuesY]',2*lendn,1)];


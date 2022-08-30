function StreamLines(X,velo,npx,npy,sx,sy,fig_num)
% StreamLines(X,velo,npx,npy,sx,sy,fig_num)
% Plot of streamlines of velocity field velo
% using Matlab's function streamline
%

if nargin < 4
    error('Too few input arguments')
elseif nargin == 4
    sx = [];
    sy = [];
    fig_num = [];
elseif nargin == 5
    sy = [];
    fig_num = [];
elseif nargin == 6
    fig_num = [];
end

if isempty(sx)
    sx = [0:0.1:1, 0.5*ones(1,11), 0:0.1:1, 1:-0.1:0];
end
if isempty(sy)
    sy = [0.5*ones(1,11), 0:0.1:1, 0:0.1:1, 0:0.1:1];
end
  
    
CoorX = zeros(npy,npx);
CoorY = zeros(npy,npx);
U = zeros(npy,npx);
V = zeros(npy,npx);

x = X(1:npx,1);
y = X(1:npx:(npx)*(npy),2);
for i=1:npy
    CoorX(i,:) = x';
end
for i=1:npx
    CoorY(:,i) = y;
end

for i=1:npy
    U(i,:) = velo((i-1)*(npx)+1: i*(npx),1)';
    V(i,:) = velo((i-1)*(npx)+1: i*(npx),2)';
end

if isempty(fig_num)
    figure('Name','Streamlines','NumberTitle','off'); clf; 
else
    figure(fig_num); clf;
end
aux = streamline(CoorX,CoorY,U,V,sx,sy);
box on; axis equal tight;
set(gca,'FontSize',12,...
    'XTick',[0:0.25:1], 'YTick',[0:0.25:1]);
axis([0,1,0,1])

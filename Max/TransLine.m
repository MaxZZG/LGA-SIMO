
close all
clear
clc

%CtrlPts = zeros(2, 4);
CtrlPts(1 : 2, 1) = [0; 1];
CtrlPts(1 : 2, 2) = [1; 2];
CtrlPts(1 : 2, 3) = [2; 1];
CtrlPts(1 : 2, 4) = [3; 4];
CtrlPts(1 : 2, 5) = [4; 2];

% knot vector
%KntVect = [0 0 0 0 1 2 3 4 4 4 4];
KntVect = [0 0 0 0 0.5 1 1 1 1];
KntVect = KntVect ./ max(KntVect);

NCtrlPts = size(CtrlPts, 2); % number of control points
p = numel(KntVect) - NCtrlPts - 1; % order of basis functions
ParaPts = linspace(0, 1, 101); % paramatric points

Idx = FindSpan(NCtrlPts, p, ParaPts, KntVect);
N0 = BasisFuns(Idx, ParaPts, p, KntVect);

C = zeros(2, numel(ParaPts));
for i = 1 : p + 1
    C = C + bsxfun(@times, N0(:, i)', CtrlPts(:, Idx - p + i - 1));
end

figure
hold on
daspect([1 1 1])
axis equal
axis on
% Plot curve
plot(C(1, :), C(2, :),'g');
% Plot control polygon
plot(CtrlPts(1, :), CtrlPts(2, :), 'k--');
% Plot control points
plot(CtrlPts(1, :), CtrlPts(2, :),...
    'r.','MarkerSize',15);

hold on 

n = 4;
m = 0;
CP = zeros(n+1,m+1,4);

CP(1,1,:) = [0,1,0,1];
CP(2,1,:) = [1,2,0,1];
CP(3,1,:) = [2,1,0,1];
CP(4,1,:) = [3,4,0,1];
CP(5,1,:) = [4,2,0,1];

Tsurface = T_Surface(CP);
Tsurface.plotSurf;

daspect([1 1 1])
axis equal

clear C CP i Idx KntVect m n N0 NCtrlPts ParaPts;
%% 使用最小二乘法拟合这条曲线

% 从T样条转换到NURBS
% 取采样点，每个节点区间取三个均匀点
numSamp = 3;
numU = (numel(Tsurface.U) - 1) * numSamp;
sampU = linspace(Tsurface.U(1+2).Value,Tsurface.U(end-2).Value,numU);
sampU = linspace(0.125,0.875,numU);
% 去掉第一个采样点和最后一个采样点
% 去掉第一个采样点和最后一个采样点
sampU = sampU(2:end-1);

% 布置NURBS的节点向量
NURBS.U = [0 0 0 0 0.5 1 1 1 1];

% 进行曲线拟合，端点插值，所以端点控制点的坐标已经确定
% 在NURBS的节点向量空间内取相同数量的采样点，并去掉第一个和最后一个
sampNU = linspace(0,1,numU);
sampNU = sampNU(2:end-1);

clear numSamp numU;

% 计算采样点处的基函数值

n = length(NURBS.U) - p - 1;
m = length(sampNU);
N = zeros(m,n - 2);

for j = 1:m
    for i = 1:n-2
        [Nip] = NURBSbasis(i, p, sampNU(j), NURBS.U, ...
                ones(length(NURBS.U) - p - 1));
        N(j,i) = Nip;
    end
end

% 计算列向量

R = zeros(m,3);

[x,y,z] = Calculate(Tsurface,0.125,0);
Q0 = [x,y,z];
[x,y,z] = Calculate(Tsurface,0.875,0);
Qend = [x,y,z];

for j = 1:m
     [x,y,z]= Calculate(Tsurface,sampU(j),0);
     plot(x,y,'.');
     hold on

     [N0] = NURBSbasis(0, p, sampNU(j), NURBS.U, ...
                ones(length(NURBS.U) - p - 1));
     [Nn] = NURBSbasis(n-1, p, sampNU(j), NURBS.U, ...
                ones(length(NURBS.U) - p - 1));
     R(j,:) = [x,y,z] - N0 * Q0 - Nn * Qend;

end

Rall = zeros(n-2,3);
for i = 1:n-2
    for j = 1:m
         [N0] = NURBSbasis(i, p, sampNU(j), NURBS.U, ...
            ones(length(NURBS.U) - p - 1));
         Rall(i,:) = Rall(i,:) + N0 * R(j,:);
    end
end


% 计算控制点：

P =  Rall / (N' * N);


clNew = [Q0(1,1:2);P(:,1:2);Qend(1,1:2)];
% clNew = sortrows(clNew)';
clNew = clNew';

KntVect = [0 0 0 0 0.5 1 1 1 1];
KntVect = KntVect ./ max(KntVect);

NCtrlPts = size(CtrlPts, 2); % number of control points
p = numel(KntVect) - NCtrlPts - 1; % order of basis functions
ParaPts = linspace(0, 1, 101); % paramatric points

Idx = FindSpan(NCtrlPts, p, ParaPts, KntVect);
N0 = BasisFuns(Idx, ParaPts, p, KntVect);

C = zeros(2, numel(ParaPts));
for i = 1 : p + 1
    C = C + bsxfun(@times, N0(:, i)', clNew(:, Idx - p + i - 1));
end

figure
hold on
daspect([1 1 1])
axis equal
axis on
% Plot curve
plot(C(1, :), C(2, :),'r');
% Plot control polygon
plot(CtrlPts(1, :), CtrlPts(2, :), 'k--');
% Plot control points
plot(CtrlPts(1, :), CtrlPts(2, :),...
    'r.','MarkerSize',15);



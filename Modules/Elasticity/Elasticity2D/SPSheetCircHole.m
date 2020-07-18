% Sheet with Small Circular Hole
%{
Copyright (C) <2014-2016>  <Khanh Chau-Nguyen, Hung Nguyen-Xuan>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}

close all
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------PRE-PROCESSING------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% control points
tic     % 开始计时
a = 1;  % 圆弧半径
L = 5;  % 四边形边长

CtrlPts1 = zeros(4, 4); 
CtrlPts1(1 : 3, 1) = [0; L; 0];
CtrlPts1(1 : 3, 2) = [-L; L; 0];
CtrlPts1(1 : 3, 3) = [-L; L; 0];
CtrlPts1(1 : 3, 4) = [-L; 0; 0];

CtrlPts1(4, :) = 1;

Curv2 = NURBSCirc((a + L) / 2, [0, 0, 0], pi / 2, pi);  % 中间的圆弧
Curv2 = HRefine(Curv2, 1, 1/2);  
CtrlPts2 = Curv2.CtrlPts4D;  
clear Curv2

Curv3 = NURBSCirc(a, [0, 0, 0], pi / 2, pi);  % 边上的圆弧
Curv3 = HRefine(Curv3, 1, 1/2);
CtrlPts3 = Curv3.CtrlPts4D;
clear Curv3

KntVect{1} = [0 0 0 1/2 1 1 1];
KntVect{2} = [0 0 0 1 1 1];

CtrlPtsSurf = cat(3, CtrlPts1, CtrlPts2, CtrlPts3);  % 连结数组，控制顶点变成(4,4,3)的三位存储结构
Surf = CreateNURBS(KntVect, CtrlPtsSurf);  

Surf = KRefine(Surf, [10, 10], [4, 4], [3, 3]);

figure
hold on
axis equal
daspect([1, 1, 1])
axis off
PlotGeo(Surf);
PlotKnts(Surf);
PlotCtrlPts(Surf);
PlotCtrlNet(Surf);
% --------------------------------------------------------------
Mesh = Mesh2D(Surf, 'VectorField');

% Assemble siffnesss matrix
disp([num2str(toc),'  Assembling the system'])
% material properties
E = 10 ^ 5;
nu = 0.3;

% 计算刚度矩阵 size(KVals) = [(NEN * 2) ^ 2, NEL]，NEN是局部基函数的数量，NEL是网格的数量
KVals = calcLocalStiffnessMatrices2D(Mesh, Surf, E, nu, 'PlaneStress');  
[Rows, Cols, Vals] = convertToTripletStorage(Mesh, KVals); % 将其转换成三元组存储方法

% Convert triplet data to sparse matrix
K = sparse(Rows, Cols, Vals);  % 这个意思是三元组中值为 0 的元素不存储了，节省空间
clear Rows Cols Vals

f = zeros(Mesh.NDof, 1);  % Mesh.NDof 是整个算例的自由度，此例中Mesh.NDof = NCtrlPts * 2

% Impose natural boundary conditions
disp([num2str(toc),'  Imposing natural boundary conditions'])

Tx = 1; 
S = SPSheetCircHoleExactStress(a, Tx);

Sxx = @(x, y) -S.xx(x, y) * (abs(x + L) <= 1e-6);  
% 边界条件是加在了一条边上，(abs(x + L) <= 1e-6)这条语句令不在 X =L的边界条件等于0
% 这样就只会对想加边界条件的 Mesh.boundry 的一部分 加边界条件了
Sxy = @(x, y) -S.xy(x, y) * (abs(x + L) <= 1e-6);
Syy = @(x, y) S.yy(x, y) * (abs(y - L) <= 1e-6);
Syx = @(x, y) S.xy(x, y) * (abs(y - L) <= 1e-6);

[Fx1, DofsFx1] = applyNewmannBdryVals(Surf, Mesh, Sxx, 3, 'FX'); % 边界依然是3号边界，但是3号边界上 X != L 部分的Sxx = 0，相当于不施加边界条件
[Fy1, DofsFy1] = applyNewmannBdryVals(Surf, Mesh, Sxy, 3, 'FY');

[Fx2, DofsFx2] = applyNewmannBdryVals(Surf, Mesh, Syx, 3, 'FX');
[Fy2, DofsFy2] = applyNewmannBdryVals(Surf, Mesh, Syy, 3, 'FY');

f(DofsFx1) = f(DofsFx1) + Fx1;
f(DofsFy1) = f(DofsFy1) + Fy1;

f(DofsFx2) = f(DofsFx2) + Fx2;
f(DofsFy2) = f(DofsFy2) + Fy2;

% Impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
h = @(x, y) 0;

[UY, DofsY] = projDrchltBdryVals(Surf, Mesh, h, 2, 'UY');
[UX, DofsX] = projDrchltBdryVals(Surf, Mesh, h, 1, 'UX');

BdryIdcs = [DofsY; DofsX];
BdryVals = [UY; UX];

FreeIdcs = setdiff(1 : Mesh.NDof, BdryIdcs);

d = zeros(Mesh.NDof, 1);
d(BdryIdcs) = BdryVals;

f(FreeIdcs) = f(FreeIdcs) - K(FreeIdcs, BdryIdcs) * BdryVals; 

% Solve the system
disp([num2str(toc),'  Solving the system'])
d(FreeIdcs) = K(FreeIdcs, FreeIdcs) \ f(FreeIdcs);

% compare stress between IGA and exact solution
% -----------------------------------------------------------------------
% IGA solution
[C, Stress] = StressEval(Surf, d, {linspace(0, 1, 401), 1}, E, nu,...
'PlaneStress');

x = reshape(C(1, :, :), 1, [])';
y = reshape(C(2, :, :), 1, [])';
alfa = atan2(y, x);

figure
hold on
plot(rad2deg(alfa), S.xx(x, y), 'r', 'LineWidth', 1.5);
plot(rad2deg(alfa), reshape(squeeze(Stress(1, :, :)), [], 1), 'b')
% -----------------------------------------------------------------------

% Export solution to vtk format
disp([num2str(toc),'  Exporting result to *.vtk format'])
ParaPts = {setdiff(linspace(0, 1, 201), 0.5), linspace(0, 1, 201)};
    
SPToVTKStress(Surf, d, ParaPts, E, nu, 'PlaneStress',...
    'SPSheetCircHoleStress')

filename = 'SPSheetCircHoleDispl';
fieldname = 'Displ';
SPToVTK(Surf, d, ParaPts, filename, fieldname)


function [f, d, FreeIdcs] = applyDrchltBdryVals(BdryIdcs, BdryVals, K, f)
% function [f, d, freeIdcs] = applyDrchltBdryVals(BdryIdcs, BdryVals, K, f)
% Apply boundary values for Dirichlet boundary condition
% 施加狄利克雷边界条件，例如位移、扭转等。在数学中狄利克雷边界条件是指微分方程的解在边界处的值。
% 
% ----------------------------------------------------------------
% Input: 
%       BdryIdcs: indices of boundary control points   要操作的边界控制点的索引
%       BdryVals: corresponding values at these control points    这些控制点对应的值，通过最小二乘拟合方法算出的值
%       K: global stiffness matrix     总刚矩阵
%       f: force vector (right hand side)  力向量（右手坐标系）
% Output:
%       f: modified force vector  修正的力向量
%       d: displacement vector    位移向量
%       freeIdcs: vector contains indices of unknowns   自由控制点索引?
% -----------------------------------------------------------------

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

NDof = numel(f);
FreeIdcs = setdiff(1 : NDof, BdryIdcs);		% 除了边界控制点，其他的控制点
d = zeros(NDof, 1);		% 位移的列向量
d(BdryIdcs) = BdryVals;		% 边界条件直接加载控制点吗
f(FreeIdcs) = f(FreeIdcs) - K(FreeIdcs, BdryIdcs) * BdryVals;  % 
end
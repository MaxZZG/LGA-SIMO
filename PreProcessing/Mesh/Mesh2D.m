function Mesh = Mesh2D(NURBS, lab)
% Mesh = Mesh2D(NURBS, type)
% modified from igafem package (Vinh Phu Nguyen)

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


if strcmp(lab, 'ScalarField') % 为什么标量域的Dof为1，以下同问，
    Dof = 1; % number of degree of freedom
elseif strcmp(lab, 'VectorField') 
    Dof = 2;
elseif strcmp(lab, 'Plate')
    Dof = 1;
else
    error('"lab" must be "ScalarField", "VectorField" or "Plate"');
end

NElDir = zeros(1, NURBS.Dim); % number of elements per direction
for i = 1 : NURBS.Dim
    NElDir(i) = numel(NURBS.uqKntVect{i}) - 1;
end

% chan =
%        1 2 3 4
%        5 6 7 8
%        9 10 11 12
%        13 14 15 16
% for a 4x2x2 control points  (4x4)

chan  = zeros(NURBS.NCtrlPts(2), NURBS.NCtrlPts(1));  % 对控制点顺序编号，沿1方向先走


count = 1;
for i = 1 : NURBS.NCtrlPts(2)
    for j = 1 : NURBS.NCtrlPts(1)
        chan(i, j) = count;
        count = count + 1;
    end
end

% element connectivity in each direction
ElDir{1} = BuildConn(NURBS.Order(1), NURBS.KntVect{1});
ElDir{2} = BuildConn(NURBS.Order(2), NURBS.KntVect{2});

% combine info from two directions to build the elements
% element is numbered as follows
%  5 | 6 | 7 | 8
% ---------------
%  1 | 2 | 3 | 4
% for a 4x2 mesh

NEl = prod(NElDir);
NEN = prod(NURBS.Order + 1); % number of local basis functions

El = zeros(NEl, NEN);

iE = 1;
for iEta = 1 : NElDir(2) %遍历2方向单元
    EtaConn = ElDir{2}(iEta, :); %单元控制点编号数组，又名（Connectivity Array p13 user_guide）
    for iXi = 1 : NElDir(1) %遍历1方向单元
        XiConn = ElDir{1}(iXi, :);%1方向单元 ixi 的控制点编号
        c = 1;
        for i = 1 : length(EtaConn)
            for j = 1 : length(XiConn)
                El(iE, c) = chan(EtaConn(i), XiConn(j)); % 映射到整体的控制点编号
                c = c + 1;
            end
        end
        iE = iE + 1;
    end
end
% 网格所包含的参数，定义了网格的形状，所以不包括控制点的坐标，控制点坐标定义在NURBS
Mesh.Dof = Dof; % 2D网格单点的自由度，

Mesh.El = El; % element connection
Mesh.NEl = NEl; % number of elements
Mesh.NEN = NEN; % number of local basis functions
Mesh.NElDir = NElDir; % number of element per direction
Mesh.NDof = NURBS.NNP * Dof; % number of dof 可以看出Dof代表每个点的自由度，空间坐标表示下，自由度最高为3

mcp = NURBS.NCtrlPts(1);
ncp = NURBS.NCtrlPts(2);

%网格边界n（1234）上的所有控制点的编号 可能不对，CompDofs 啥意思？

Mesh.Boundary(1).CompDofs{1} = sub2ind([mcp, ncp], ones(1, ncp), 1 : ncp)';
Mesh.Boundary(2).CompDofs{1} = sub2ind([mcp, ncp], mcp * ones(1, ncp), 1 : ncp)';
Mesh.Boundary(3).CompDofs{1} = sub2ind([mcp, ncp], 1 : mcp, ones(1, mcp))';
Mesh.Boundary(4).CompDofs{1} = sub2ind([mcp, ncp], 1 : mcp, ncp * ones(1, mcp))';

%网格边界n（1234）下一层上的所有控制点的编号

Mesh.Boundary(1).NextLayerDofs.CompDofs{1} = sub2ind([mcp, ncp], 2 * ones(1, ncp), 1 : ncp)';
Mesh.Boundary(2).NextLayerDofs.CompDofs{1} = sub2ind([mcp, ncp], (mcp - 1) * ones(1, ncp), 1 : ncp)';
Mesh.Boundary(3).NextLayerDofs.CompDofs{1} = sub2ind([mcp, ncp], 1 : mcp, 2 * ones(1, mcp))';
Mesh.Boundary(4).NextLayerDofs.CompDofs{1} = sub2ind([mcp, ncp], 1 : mcp, (ncp - 1) * ones(1, mcp))';

for iSide = 1 : 4
    ind = mod(floor((iSide + 1) / 2), 2) + 1; % 判断NURBS参数方向
    Mesh.Boundary(iSide).NDof = NURBS.NCtrlPts(ind) * Dof;
    
    if Dof == 2 % vector field % 以下不太懂，问题回到Dof = 2 ？
        % 有两个自由度分量，前一个自由度分量排在前面，后一个自由度分量排在后面
        Mesh.Boundary(iSide).Dofs = [Mesh.Boundary(iSide).CompDofs{1}; Mesh.Boundary(iSide).CompDofs{1} + NURBS.NNP];
        Mesh.Boundary(iSide).CompDofs{2} = Mesh.Boundary(iSide).CompDofs{1} + NURBS.NNP;
        Mesh.Boundary(iSide).NextLayerDofs.CompDofs{2} = Mesh.Boundary(iSide).NextLayerDofs.CompDofs{1} + NURBS.NNP;
    elseif Dof == 1 % scalar field or plate problem
        Mesh.Boundary(iSide).Dofs = Mesh.Boundary(iSide).CompDofs{1};
    end
    Mesh.Boundary(iSide).Dof = Dof;
    Mesh.Boundary(iSide).NEl = NElDir(ind);
    Mesh.Boundary(iSide).NElDir = NElDir(ind);
    Mesh.Boundary(iSide).NEN = NURBS.Order(ind) + 1;
    Mesh.Boundary(iSide).El = ElDir{ind}; 
end

end





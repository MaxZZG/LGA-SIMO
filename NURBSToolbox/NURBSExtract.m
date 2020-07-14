function ONURBS = NURBSExtract(INURBS, Dir, val)
% function ONURBS = NURBSExtract(INURBS, Dir, val)
% ------------------------------------------------------------------
% Extract lower dimensional NURBS object.
% 提取低维度的NURBS对象?
%-------------------------------------------------------------------
% Input:
%       INURBS: input NURBS structure
%       Dir: direction along which to extract
%       (xi: dir = 1, eta: dir = 2, zeta: dir = 3)
%       val: parametric value from which to exact
%-------------------------------------------------------------------
% Output:
%       ONURBS: output NURBS structure
%-------------------------------------------------------------------

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

% if 'Dir'th dimension of CtrlPts4D matrix is not the last
% dimension, then permute it until it lies at last dimension



if Dir ~= INURBS.Dim  % 假设 Dir = 1，Dim = 3
    Dirs = 1 : INURBS.Dim + 1;  % Dirs = [1,2,3,4];
    Dirs(Dir + 1) = []; %Dirs = [1,3,4];
    Dirs = [Dirs, Dir + 1]; %Dirs = [1,3,4,2];
    temp = permute(INURBS.CtrlPts4D, Dirs);  % 对N维数组重新排列其维数，
else
    temp = INURBS.CtrlPts4D;
end
dim = size(temp);
temp = reshape(temp, [], dim(end));

CtrlPts = CurvPntByCornerCut(INURBS.NCtrlPts(Dir),...
    INURBS.Order(Dir), INURBS.KntVect{Dir}, temp, val);  %通过割角方式计算控制点,

CtrlPts = reshape(CtrlPts, dim(1 : end - 1));

KntVect = INURBS.KntVect;
KntVect{Dir} = [];
KntVect = KntVect(~cellfun(@isempty, KntVect));
ONURBS = CreateNURBS(KntVect, CtrlPts);
end
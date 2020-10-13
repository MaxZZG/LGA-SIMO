function ONURBS = KRefine(INURBS, NEl, Order, Continuity)
% function ONURBS = KRefine(INURBS, NEl, Order, Continuity)
% -----------------------------------------------------------------
% Refine a NURBS object (curve, surface, volume) by knot refinement
% and order (degree) elevation
% -----------------------------------------------------------------
% Input:
%       INURBS: input NURBS structure
%       NEl: number of elements per direction % 每个单元再细化为多少个单元
%       Order: order of basis functions per direction
%       Continuity: continuity of basis functions per direction
% -------------------------------------------------------------
% Output:
%       ONURBS: output NURBS structure
% -------------------------------------------------------------

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

ONURBS = INURBS;
% order elevation
assert(all(Order - ONURBS.Order >= 0));
if all(Order)
    for Dir = 1 : ONURBS.Dim
        p = Order(Dir);
        t = p - ONURBS.Order(Dir);
        if t > 0
            ONURBS = PRefine(ONURBS, Dir, t);
        end
    end
end
% knot refinement
assert(all(Continuity >= 0));
assert(all(Continuity < Order));
assert(all(NEl > 0));
for Dir = 1 : ONURBS.Dim
    Knts = ONURBS.uqKntVect{Dir};
    Mlts = ONURBS.KntMult{Dir};
    N = NEl(Dir);
    % compute the number of knots to be inserted
    deltaXi = diff(Knts) / N;% diff相邻元素的差分
    Knts = Knts(1 : end - 1);
    
    Xi = repmat(Knts', 1, N);% N个knts列向量
    step = 1 : N - 1;
    Xi(:, 2 : end) = Xi(:, 2 : end) + deltaXi' * step;%每一行，为每个单元待插入的节点
    
    Mlts = Mlts(1 : end - 1);
    repsMat = zeros(numel(Mlts), N); 
    repsMat(:, 1) = ONURBS.Order(Dir) - Continuity(Dir) - Mlts;% 原节点的重复度
    repsMat(:, 2 : end) = ONURBS.Order(Dir) - Continuity(Dir); % 新插入的节点的重复度
    
    Knts = reshape(Xi', 1, []);
    repsVect = reshape(repsMat', 1, []);
    
    Knts = Knts(2 : end);
    repsVect = repsVect(2 : end);
    
    insKnts = Knts(repsVect > 0); % inserted knots
    repsVect = repsVect(repsVect > 0); % 构造出新的非重节点向量，和重复度向量
    
    Idx = zeros(1, sum(repsVect));
    Idx(cumsum([1 repsVect(1 : end - 1)])) = 1;% 设置每个节点的起点
    Idx = cumsum(Idx);% 构造索引
    
    if ~isempty(insKnts)
        KntsMult = insKnts(Idx); % 获取含重的节点向量
        % insert multiple knots
        ONURBS = HRefine(ONURBS, Dir, KntsMult);
    end
end
end

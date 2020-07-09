function B = AllBernstein(p, xi)
% function B = AllBernstein(p, xi)
%--------------------------------------------------------------------------
% Compute all pth-order Bernstein basis functions.
% Bezier曲线的基函数：Bernstein 多项式
%--------------------------------------------------------------
% Input:
%      p: order of polynominal  
%      xi: parametric points  xi是一个数组，意思是多个参数u
%--------------------------------------------------------------
% Output:
%      B: bernstein basis functions B为矩阵，每个参数u对应 p+1 个Bernstein的值
%--------------------------------------------------------------
% Based on Algorithm A1.3 [The NURBS BOOK, p.20]
%--------------------------------------------------------------------------

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

n = p + 1;
B = zeros(numel(xi), n);
Bi = zeros(n, 1); 	% Bi是行向量或者列向量并不影响，最后赋值给B了。感觉Bi最好是行向量。
for jj = 1 : numel(xi)	% 参数u的个数为循环次数
    Bi(1) = 1;
    u1 = 1 - xi(jj);
    for j = 2 : n
        saved = 0;
        for k = 1 : j - 1
            temp = Bi(k);
            Bi(k) = saved + u1 * temp;
            saved = xi(jj) * temp;
        end
        Bi(j) = saved;
    end
    B (jj, :) = Bi; % B 第jj行 = Bi Bi为列向量，赋值时会将一整列赋值给B规定的行里。
end
end

function B = AllBernstein(p, xi)
% function B = AllBernstein(p, xi)
%--------------------------------------------------------------------------
% Compute all pth-order Bernstein basis functions.
% Bezier���ߵĻ�������Bernstein ����ʽ
%--------------------------------------------------------------
% Input:
%      p: order of polynominal  
%      xi: parametric points  xi��һ�����飬��˼�Ƕ������u
%--------------------------------------------------------------
% Output:
%      B: bernstein basis functions BΪ����ÿ������u��Ӧ p+1 ��Bernstein��ֵ
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
Bi = zeros(n, 1); 	% Bi����������������������Ӱ�죬���ֵ��B�ˡ��о�Bi�������������
for jj = 1 : numel(xi)	% ����u�ĸ���Ϊѭ������
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
    B (jj, :) = Bi; % B ��jj�� = Bi BiΪ����������ֵʱ�Ὣһ���и�ֵ��B�涨�����
end
end

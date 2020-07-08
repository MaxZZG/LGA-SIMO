function [AeqIn] = incompresibility(sdof,Bstrain,totalGP)
% function [AeqIn] = incompresibility(sdof,Bstrain,totalGP)
%{
Copyright (C) <2014-2016>  <Hung Nguyen-Xuan, Khanh Chau-Nguyen>

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

% incompresibility condition at every gauss point
AeqIn = sparse(totalGP,sdof);
for i = 1:totalGP
    Bx = Bstrain{i}(1,:); 
    By = Bstrain{i}(2,:);
    AeqIn(i,:) = Bx + By; 
   clear Bx By;
end

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

clear
close all
clc

radius = 2;
CtrlPts = zeros(4, 9, 2);

CtrlPts(1 : 3, 1, 2) = [1; 0; 0];
CtrlPts(1 : 3, 2, 2) = [1; 1; 0];
CtrlPts(1 : 3, 3, 2) = [0; 1; 0];
CtrlPts(1 : 3, 4, 2) = [-1; 1; 0];
CtrlPts(1 : 3, 5, 2) = [-1; 0; 0];
CtrlPts(1 : 3, 6, 2) = [-1; -1; 0];
CtrlPts(1 : 3, 7, 2) = [0; -1; 0];
CtrlPts(1 : 3, 8, 2) = [1; -1; 0];
CtrlPts(1 : 3, 9, 2) = [1; 0; 0];

CtrlPts(1 : 3, :, :) = CtrlPts(1 : 3, :, :) * radius;

% assign the weights
CtrlPts(4, :, :) = 1;
W = 1 ;
%W = 1 / sqrt(2);

CtrlPts(:, 2 : 2 : 8, 2) = CtrlPts(:, 2 : 2 : 8, 2) * W;

% knot vector
KntVect{1} = [0 0 0 1 1 2 2 3 3 4 4 4];
KntVect{1} = KntVect{1} ./ max(KntVect{1});
KntVect{2} = [0 0 1 1];

Surf = CreateNURBS(KntVect, CtrlPts);

figure
hold on
grid on
set(gcf,'color','white')
daspect([1 1 1])
axis([-radius - 0.2 radius + 0.2 -radius - 0.2 radius + 0.2]);
axis equal
axis off
PlotGeo(Surf)
PlotKnts(Surf)
PlotCtrlPts(Surf)
PlotCtrlNet(Surf)

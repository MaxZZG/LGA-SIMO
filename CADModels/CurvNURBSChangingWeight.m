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

% control points

CtrlPts = zeros(5, 4);
CtrlPts(1 : 2, 1) = [0; 1];
CtrlPts(1 : 2, 2) = [1; 2];
CtrlPts(1 : 2, 3) = [2; 1];
CtrlPts(1 : 2, 4) = [3; 4];
CtrlPts(1 : 2, 5) = [4; 2];

CtrlPts(4, :) = 1;

W = [1, 1, 1, 1];
Knts = [0 0 0 0 0.5 1 1 1 1];

Curv = cell(numel(W), 1);

for i = 1 : numel(W)
    temp = CtrlPts;
    temp(:, 2) = temp(:, 2) .* W(i);
    Curv{i} = CreateNURBS({Knts}, temp);
end

figure
hold on
axis off
set(gcf,'color','white')
set(gca,'XTick', 0 : 1: 3)
set(gca,'YTick', 0 : 1: 2)

for i = 1 : numel(W)
    PlotGeo(Curv{i})
    PlotCtrlPts(Curv{i})
    PlotKnts(Curv{i})
    PlotCtrlNet(Curv{i})
end


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

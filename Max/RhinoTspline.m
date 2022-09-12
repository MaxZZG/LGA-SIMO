

clc
clear all
close all

Tspline = read_bezier_extraction('plane.iga');

drawConPts(Tspline.control_points);

% 计算出所有的Bezier单元


noElem = numel(Tspline.elements);
Q = zeros(16,3,noElem);
for e=1:noElem
    el           = Tspline.elements(e);
    element{e}   = el.connectivity;
    C{e}         = el.extraction;
    degrees(e,:) = el.degree;


    we     = diag(ones(16,1));% element weights
    Wb     = C{e}'* ones(16,1);  % element Bezier weights
    Q(:,:,e) = inv(diag(Wb)) * C{e}' * we * Tspline.control_points(1:3,element{e})';
end

% 生成NURBS节点向量

uknot = [0,0,0,0,1.5,2,2,2,2];
vknot = [0,0,0,0,1.5,2,2,2,2];

uknot = uknot / max(uknot);
vknot = uknot / max(vknot);
p = 3;
q = 3;

C_NURBS = bezierExtraction2D(uknot,vknot,p,q);

% 生成patchCAE数据结构

noElems = 4;
% 根据Bezier曲面计算NURBS的控制点
elements = zeros(4,16);
elements(1,:) = [1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19];
elements(2,:) = [2,3,4,5,7,8,9,10,12,13,14,15,17,18,19,20];
elements(3,:) = [6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24];
elements(4,:) = [7,8,9,10,12,13,14,15,17,18,19,20,22,23,24,25];

for e=1:noElems
    
    Ce     = C_NURBS(:,:,e);             % element Bezier extraction operator
    we     = diag(ones(16,1));% element weights
    Wb     = Ce'*ones(16,1);  % element Bezier weights
    
    bezierPts = Q(:,:,e);

    A = inv(diag(Wb))*Ce'*we;
    pts = inv(A) * bezierPts;

    test = inv(diag(Wb))*Ce'*we * pts;

end




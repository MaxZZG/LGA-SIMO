
clc;
clear all;
close all;



% NURBS bezier-extraction

Nknot = [0,0,0,0,1,2,3,4,4,4,4];
[B,nb] = bezierExtraction(Nknot,3);

Nknot = [-2,-2,-2,-2,-1,0,1,2,2,2,2];
[B1,nb] = bezierExtraction(Nknot,3);

knot = 0:4;
U = [1,1,1,2,2,2,3,3,3];
spans = [1,1,1,2,2,2,3,3,3];

U = [1,2,3];
spans = [1,2,3];
m = numel(U);
A = TsplineExtraction(knot,U,spans,3,m);




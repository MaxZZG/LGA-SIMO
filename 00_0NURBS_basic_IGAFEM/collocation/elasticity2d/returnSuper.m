function [ sy, sdy, sddy ] = returnSuper( p )
%returns the superconvergent points in the reference interval [-1,1]

%define some constants
sp0 = 1/sqrt(3);
sp1 = (1/15)*sqrt(225-30*sqrt(30));
sp2 = 0.50491856751265330603;
sp3 = 0.503221894597504;

switch p
    case 1
        sy = [-1,1];
        sdy = 0;
        sddy = [];
    case 2
        sy = [-1, 0, 1];
        sdy = [-sp0, sp0];
        sddy = 0;
    case 3
        sy = [-sp1, sp1];
        sdy = [-1, 0, 1];
        sddy = [-sp0, sp0];
    case 4
        sy = [-1, 0, 1];
        sdy = [-sp1, sp1];
        sddy = [-1, 0, 1];
    case 5
        sy = [-sp2, sp2];
        sdy = [-1, 0, 1];
        sddy = [-sp1, sp1];
    case 6
        sy = [-1, 0, 1];
        sdy = [-sp2, sp2];
        sddy = [-1, 0, 1];
    case 7
        sy = [-sp3, sp3];
        sdy = [-1, 0, 1];
        sddy = [-sp2, sp2];
    otherwise
        disp('Invalid p')
        sy = [];
        sdy = [];
        sddy = [];                
end


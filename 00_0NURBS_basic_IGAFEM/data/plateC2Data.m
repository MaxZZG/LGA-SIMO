% rectangular plate in tension



controlPts = [
0   0;
0.3 0;
0.6 0;
1.0 0;
0.   0.3;
0.3  0.3;
0.6  0.3;
1    0.3;
0    0.6;
0.3  0.6;
0.6  0.6;
1.0  0.6;
0 1;
0.3 1;
0.6 1;
1 1];

uKnot = [0 0 0 0 1 1 1 1];
vKnot = [0 0 0 0 1 1 1 1];

p = 3;
q = 3;

noPtsX = 4;
noPtsY = 4;

weights = ones(1,noPtsX*noPtsY)';


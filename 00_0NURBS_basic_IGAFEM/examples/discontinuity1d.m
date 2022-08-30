addpath('../C_files')
%addpath('~/code/xfem-efg-matlab/fem_util');

clear all

% original data

knots      = [0 0 0 1 1 1];
controlPts = [0 0; 0.5 0.5; 1 0];
p          = 2;
noPtsX     = 3;
weights    = [1 1 1]; 

% plot basis and curve

noPts       = 60;
xi          = linspace(0,1,noPts);
BsplineVals = zeros(noPts,noPtsX);

for i=1:noPtsX
    for c=1:noPts
        [BsplineVals(c,i) NURBSderivs(c,i+1)] = ...
            NURBSbasis(i, p, xi(c), knots, weights);
    end
end

sCurve      = zeros(2,noPts);
nCurve      = zeros(2,noPts);

for i=1:noPts
    sCurve(1,i) = NURBSinterpolation(xi(i), p, knots, ...
                                controlPts(:,1), weights);
    sCurve(2,i) = NURBSinterpolation(xi(i), p, knots, ...
                                controlPts(:,2), weights);
end

figure(1)
plot(xi, BsplineVals,'LineWidth',1.5);

figure(2)
hold on
plot(controlPts(:,1),controlPts(:,2),'r-o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',10);
plot(sCurve(1,:),sCurve(2,:),'r-','LineWidth',1.8);

% knot insertion p+1 times to create discontinuity

newKnotsX = [1/2 1/2 1/2];

uniqueKnots = unique(knots);
nonewkX     = size(newKnotsX,2);
weightedPts = [controlPts weights'];

[newKnots,newControlPts] = ...
    RefineKnotVectCurve(noPtsX-1,p,knots,weightedPts,newKnotsX,nonewkX-1);

knots       = newKnots;
controlPts  = [newControlPts(:,1)./newControlPts(:,3) ...
               newControlPts(:,2)./newControlPts(:,3) newControlPts(:,3)];

noPtsX = noPtsX + nonewkX;

% plot the new basis

noPts1 = 150;
xi          = linspace(0,1,noPts1);
BsplineVals = zeros(noPts1,noPtsX);

for i=1:noPtsX
    for c=1:noPts1
        [BsplineVals(c,i) NURBSderivs(c,i+1)] = ...
            NURBSbasis(i, p, xi(c), knots, controlPts(:,3));
    end
end

plot(xi, BsplineVals,'LineWidth',1.5);

% plot the curve

controlPts(4,2) = controlPts(4,2)+0.01;
for i=1:noPts1
    sCurve(1,i) = NURBSinterpolation(xi(i), p, knots, ...
                                controlPts(:,1), controlPts(:,3));
    sCurve(2,i) = NURBSinterpolation(xi(i), p, knots, ...
                                controlPts(:,2), controlPts(:,3));
end

figure(3)
hold on
plot(controlPts(:,1),controlPts(:,2),'r-o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',10);
plot(sCurve(1,1:75),sCurve(2,1:75),'r-','LineWidth',1.8);
plot(sCurve(1,76:150),sCurve(2,76:150),'b-','LineWidth',1.8);
axis tight

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
exportfig(gcf,'curve1d.eps',opts)




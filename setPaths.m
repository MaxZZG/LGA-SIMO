clc
clear
close all

% NURBS basis
path('./00_0NURBS_basic_IGAFEM/fem_util', path);
path('./00_0NURBS_basic_IGAFEM/C_files_win', path);
path('./00_0NURBS_basic_IGAFEM/data', path);
path('./00_0NURBS_basic_IGAFEM/meshing', path);
path('./00_0NURBS_basic_IGAFEM/post-processing', path);
path('./00_0NURBS_basic_IGAFEM/fem-functions', path);
path('./00_0NURBS_basic_IGAFEM/analytical-solutions', path);
path('./00_0NURBS_basic_IGAFEM/nurbs-util', path);
path('./00_0NURBS_basic_IGAFEM/bezier-extraction', path);

% SimoPackage
% Add the following to the path
addpath(genpath(fullfile(pwd,'CADModels')))
addpath(genpath(fullfile(pwd,'Modules')))
addpath(genpath(fullfile(pwd,'NURBSToolbox')))
addpath(genpath(fullfile(pwd,'PostProcessing')))
addpath(genpath(fullfile(pwd,'PreProcessing')))
addpath(genpath(fullfile(pwd,'Processing')))
addpath(genpath(fullfile(pwd,'Tspline')))

addpath(genpath(fullfile(pwd,'Max'))) % Max的例子

addpath(genpath(fullfile(pwd,'input'))) % Max的输入



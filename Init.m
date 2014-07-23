% Initializes absolute paths
addpath .
addpath(genpath(fullfile(pwd,'Codes1D')))
% addpath(genpath(fullfile(pwd,'Codes2D')));rmpath(genpath(fullfile(pwd,'Codes2DQuad')));
addpath(genpath(fullfile(pwd,'Codes2DQuad')));rmpath(genpath(fullfile(pwd,'Codes2D')));
addpath(genpath(fullfile(pwd,'CodesDPG')))
addpath(genpath(fullfile(pwd,'ServiceRoutines')))
addpath(genpath(fullfile(pwd,'Grid'))) % recursive path
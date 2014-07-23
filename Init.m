% Initializes absolute paths
function Init(useQuads)
if nargin<1
    useQuads = 0; % default to triangles
end
addpath .
addpath(genpath(fullfile(pwd,'Codes1D')))
if useQuads
    addpath(genpath(fullfile(pwd,'Codes2DQuad')));rmpath(genpath(fullfile(pwd,'Codes2D'))); % comment out if you want triangles instead
else
    addpath(genpath(fullfile(pwd,'Codes2D')));rmpath(genpath(fullfile(pwd,'Codes2DQuad')));
end
addpath(genpath(fullfile(pwd,'CodesDPG')))
addpath(genpath(fullfile(pwd,'ServiceRoutines')))
addpath(genpath(fullfile(pwd,'Grid'))) % recursive path
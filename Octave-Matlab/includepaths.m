%INCLUDEPATHS Script to add the folders of each method
%   This script adds to the workspace the folders of each optimization
%   method, the test problems and the quality indicators

cp = pwd; %gets the current path
addpath([cp '/Algorithms']); %path for the algorithms
addpath([cp '/DTLZ']); %path for the DTLZ functions

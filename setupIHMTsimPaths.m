function baseDir = setupIHMTsimPaths

% Edit this function to contain the paths of your folders:

%% SimWithSpoil Code base directory:
baseDir = '/path/to/OptimizeIHMTimaging/';
addpath(genpath(baseDir)) 
cd baseDir; % used for easy saving later

%% qMRlab code directory:
% http://qmrlab.org/
% https://github.com/qMRLab/qMRLab

qMRlabDirectory = '/path/to/qMRLab-master';
addpath(genpath(qMRlabDirectory)) 

%% MP2RAGE code:
% https://github.com/JosePMarques/MP2RAGE-related-scripts
mp2rageDirectory = '/path/to/MP2RAGE-related-scripts';
addpath(genpath(mp2rageDirectory)) 

%% Misc code:
% https://github.com/christopherrowley/NeuroImagingMatlab
neuroMatlabDir = '/path/to/NeuroImagingMatlab';
addpath(genpath(neuroMatlabDir)) 




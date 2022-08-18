%% Run Test Cases:

% this script runs the different test cases, and outputs the figures in the
% according directory.

savDir =  '/Path/To/OptimizeIHMTimaging/TestCases/Figures/';

FLASH_test(savDir)

% Here we export a table of M0 and T1 values to compare
VFA_resultsTable = VFA_test(savDir);

% MPRAGE test, not sure why, just for fun
MPRAGE_test(savDir)

% With the excitation part of the simulations confirmed to work (?)...
% look to saturation part with a qMT style approach.
qMT_resultsTable = qMT_test(savDir);



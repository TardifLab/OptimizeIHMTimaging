
%% Run one more set of simulations with the fit parameters to confirm the fit. 
TParams = [ 50;...      % R
            50e-3;...   % T2a
            0.75e-3;...    % T1D
            11.5e-6;...  % T2b
            0.071;...   % M0b
           0.25];         % R1b
    
idx = 1; % used as a suffix to save the fit figure.

baseDir = '/path/to/OptimizeIHMTimaging/';
saveImgDir = '/Directory/To/Save/Output/Figures/';


SingleFitTestingFunction(baseDir, saveImgDir, TParams, idx);




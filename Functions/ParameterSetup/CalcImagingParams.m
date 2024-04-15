function Params = CalcImagingParams(Params)
%% These should be hardware related parameters that you likely wont be changing.
% Anything that can change that you might end up looping through should go
% into CalcVariableImagingParms.m

% The goal of this function is to fill the Params structure up with the
% required values for calculations. 

% Also to check and put in some default values in case values are either
% not set, and not known. 

if ~isfield(Params,'CalcVector') % turn off if you don't care about longitudinal change over time in magnetization
    Params.CalcVector = 0; % save resources
end

if ~isfield(Params,'B0') % if not defined, assume 3T
    Params.B0 = 3; % main field strength (in Tesla)
end

if ~isfield(Params,'RFspoiling')
    Params.RFspoiling = true; % default
end

if ~isfield(Params,'RFspoilingPhaseInc') % if not defined,assume spoiling in readout
    % In general, value is based on vendor GE = 117, Siemens = 50, Philips = 150 
    Params.RFspoilingPhaseInc = 50;  % in degrees 
end

if ~isfield(Params,'PerfectSpoiling') % if not defined,assume spoiling in readout
    Params.PerfectSpoiling = true; % eliminate transverse, used for testing
end

if ~isfield(Params,'GradientSpoiling') % if not defined,assume spoiling in readout
    Params.GradientSpoiling = true; % binary
    Params.PerfectSpoiling = false; % eliminate transverse, used for testing
end


if ~isfield(Params,'GradientSpoilingStrength') % if not defined,assume spoiling in readout
    Params.GradientSpoilingStrength = 24; % mT/m
end

if ~isfield(Params,'IncreasedGradSpoil') % if not defined,assume spoiling in readout
    Params.IncreasedGradSpoil = true; % binary
end

if ~isfield(Params,'Params.ModelSpinDiffusion') % if not defined, don't account for spin diffusion
    Params.ModelSpinDiffusion = true; % binary
end


if ~isfield(Params,'ReadoutResolution') % if not defined,assume spoiling in readout
    Params.ReadoutResolution = 1e-3; % m - for calculating spoiling
end

if ~isfield(Params,'A_g') || ~isfield(Params,'maxDephase') || ~isfield(Params,'G_t')
    % This is currently going to override any this set to the above three
    % paramters. You can either modify code, or back calculate a gradient
    % spoiling strength to give you the same gradient area for a set time.
    % Calculate Spoiling moment from gradient
    Params = CalculateGradientSpoilingMoment( Params, 0 );
end

if ~isfield(Params,'stepSize') % if not defined,assume spoiling in readout
    Params.stepSize = 50e-6; % m - for calculating spoiling
end

if Params.B0 == 7
    if ~isfield(Params,'TransmitCoil') % if not defined, assume 3T
        Params.TransmitCoil = 'STX'; % main field strength (in Tesla)
    end
end



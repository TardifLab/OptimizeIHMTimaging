function Params = CalcVariableImagingParams(Params)

%% Variables that might be adjusted between sequences. Set some defaults


%% Excitation Parameters
if Params.numExcitation == 1 % for standard GRE, make sure you have this.
    Params.DummyEcho = 0;
end

if ~isfield(Params,'echoSpacing') % for standard GRE, make sure you have this.
    Params.echoSpacing = 5e-3;
end

if ~isfield(Params,'WExcDur') % if not defined,assume spoiling in readout
    Params.WExcDur = 0.1/1000;  %in seconds
end

% Params.ExcB1 = Params.flipAngle /(360*42.577478518*Params.WExcDur);


%% Tissue Parameters


if ~isfield(Params,'N_spin') % if not defined,assume spoiling in readout
    Params.N_spin = 201; % mm - for calculating spoiling

    % Value independently checked to be good, but taken from: 
    % Yarnykh, V.L., 2010. Optimal radiofrequency and gradient spoiling for improved 
    % accuracy of T1 and B1 measurements using fast steady-state techniques. Magn. 
    % Reson. Med. 63, 1610–1626. https://doi.org/10.1002/mrm.22394
end

%% Saturation Parameters
if ~isfield(Params,'boosted') % if not defined, assume non-boosted protocol
    Params.boosted = 0; % binary
end

% Calculate spoiling moment for end of MTsat train. 
if ~isfield(Params,'MTC') % if not defined, no MT
    Params.MTC = 0; % binary
end

if Params.MTC
    % Allow different spoiling for MT
    if Params.GradientSpoiling 
        if ~isfield(Params,'GradientSpoilingStrength_MT')
            Params.GradientSpoilingStrength_MT = 20; % mT/m
        end
    end

    if ~isfield(Params,'A_g') || ~isfield(Params,'maxDephase') || ~isfield(Params,'G_t')
        % Calculate Spoiling moment from gradient
        Params = CalculateGradientSpoilingMoment( Params, 1 );
    end
end


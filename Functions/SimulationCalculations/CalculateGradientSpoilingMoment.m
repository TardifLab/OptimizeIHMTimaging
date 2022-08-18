function Params = CalculateGradientSpoilingMoment( Params, MTC )

% Export is Params.A_grad,  -> magnitude and phase to calculate rotation.

%% You will need to know the values for your particular scanner.

% We need the following parameters to be set:
% Params.GradientSpoilingStrength = 24; % mT/m
% Params.ReadoutResolution = 1; % mm - for calculating spoiling
% Params.IncreasedGradSpoil = true; % binary

% Find attributes related to the gradient spoiling

% Position within Voxel
spin_dis = Params.ReadoutResolution;

if ~MTC
    GradientSpoilingStrength = Params.GradientSpoilingStrength;

    if Params.IncreasedGradSpoil && ~isfield(Params,'G_t')
        Params.G_t = 2.625e-3;
    elseif ~isfield(Params,'G_t')
        Params.G_t = 1.645e-3; % in seconds
    end

else
    GradientSpoilingStrength = Params.GradientSpoilingStrength_MT;
    Params.G_t_MT = 1.1e-3; % in seconds
    Params.G_time_elapse_MT = 1.6e-3;

end

% Field experienced is mT/(m), get field by multiplying gradient by
% position
% spin_field = GradientSpoilingStrength * spin_dis; % in T  (T/m *m)

% Larmor precession:
% w = 2*pi *42.577478518E6 * spin_field; % rad/s as a function of position along gradient

% Separate to save to different variable names:
Params.A_g = GradientSpoilingStrength/1000 * Params.G_t; % Dephase Gradient area over time

if MTC
    Params.A_g_MT = GradientSpoilingStrength/1000 * Params.G_t_MT; % Dephase Gradient area over time
end



% Dephasing moment = w * G_t;













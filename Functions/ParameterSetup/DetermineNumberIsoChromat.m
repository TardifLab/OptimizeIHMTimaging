function N_spin = DetermineNumberIsoChromat(Params, TD)

%% Need to determine a sufficient number of isochromats. 
% From Gudbjartsson, H., Patz, S., 1995. IEEE Trans. Med. Imaging, 
% we need x spacing small enough such that we see some spin diffusion over
% the smallest time step. Rearrange their equation embeded in text on pg
% 639 to get:

idx = 1;
t_values = zeros(1,10);
if Params.MTC
    t_values(idx) = Params.pulseGapDur;
    idx = idx+1;
    t_values(idx) = Params.G_t_MT;
    idx = idx+1;
    if Params.boosted
        t_values(idx) = Params.TD_MT;
        idx = idx+1;
    end
end
t_values(idx) = Params.echoSpacing;
idx = idx+1;
t_values(idx) = TD;

% At 1ms or shorter, we exceed 1000 spins, at that point, we will not model diffusion 
t_values(t_values < 3e-3) = [];
t_min = min(t_values);

N_spin =ceil(sqrt(Params.ReadoutResolution^2/(2*Params.D*t_min)));

% For consistency with the diffusion convolution, lets keep number odd
N_spin = N_spin-mod(N_spin,2)+1;

% Ensure minimum number of 180 based on Varynkh 2010
if N_spin < 180
    N_spin = 180;
end
function Params = CalcBoundSatFromInversionPulse(Params)

% Reduce computational load by precomputing the impact of the excitation
% pulses on the bound pool.

% Assuming instanteous saturation, ignore exchange and relaxation.
delta = 0;

b1 = (Params.InversionEfficiency*180) /(360*42.577478518*Params.InvPulseDur);
w1 = 2*pi *42.577478518 * b1;  % assume B1 is in microTesla, and drop the 10^6 from gamma. w1 in rad/s

%% For consistency, use computeG code from qMRLab

G = computeG(delta, Params.T2b, Params.lineshape);

if strcmp(Params.lineshape, 'SuperLorentzian') % default is gaussian
    Params.Rrfb_inv = pi.*w1.^2.* G;    
elseif strcmp(Params.lineshape, 'Lorentzian') % default is gaussian
    Params.Rrfb_inv = pi.*w1.^2.* G;    
else
    Params.Rrfb_inv = w1.^2 .* G;
end


% For dipolar pool == 0 since delta = 0.
Params.Rrfd_inv = 0; % Params.Rrfb_exc*(2*pi*delta/wloc)^2;



%% There is the option to piggy back off of qMRlab code
% but this seems to work better for super short pulses.
% if you wanted to implement it though...

% Rpulse = GetPulse( flipAngle, 0, Params.WExcDur, 'sinc');
% % add all since we do it in one time step
% w1  = sum(Rpulse.omega(0:Params.stepSize:Params.WExcDur)); 

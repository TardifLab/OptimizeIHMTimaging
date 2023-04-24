function Params = CR_SAR_scale_PulseHeight(Params, B1peak_limit)
% Uses SAR restrictions to give an estimate of the maximum pulse height 
% that can be used for the saturation pulses
% necessary parameters:
% numSat = number of sat pulses
% satRMS = root mean square of sat pulses (in Tesla)
% tp_sat = time of sat pulse (in seconds)
% numExc = number of excitation pulses
% flip = flip angle of excitation pulses (in degrees) 
% TR = time (in seconds)
% B1peak_limit in microTesla


if ~isfield(Params,'B0') % if not defined, assume 3T
    Params.B0 = 3; % main field strength (in Tesla)
end

if ~isfield(Params, 'boosted')
    error ('Please specify if boosted sat scheme is used (enter Params.boosted = 0 or 1)')
end

if Params.B0 == 7
    if strcmp(Params.TransmitCoil, 'STX')
        SAR_limit = 3; % Empirical value to match what I get at scanner
    end

    empFact = -1.35e-5*Params.TR + 1.21e-3; % Rough estimate

else
    SAR_limit = 3; %(W/kg)
    empFact = 1.44e-3; % 3T
end
    
gam = 42.576e6;
kg = 60; % reference weight
w0 = gam *Params.B0;


%% Power = J/s. Multiply by pulse time to find J of work done
% For Excitation pulse
% (empiricalFactor*B1field^2) * numberPulses * time pulse
excB1 = (Params.flipAngle) / (360* gam * Params.WExcDur);
J_exc = (empFact * excB1^2 * w0^2)* (Params.numExcitation * Params.WExcDur); % Power (J/s) * time (s) = J

Jsat = SAR_limit*Params.TR*kg - J_exc; % J/(s*kg) *s*kg = J - J = J

if Jsat < 0
    Params.satRMS = 0;
else

    if Params.boosted % modify for different definition of numSatPulse
        SatPulseNumberTotal = Params.numSatPulse * Params.satTrainPerBoost;
    else
        SatPulseNumberTotal = Params.numSatPulse;
    end
    
    % Power sat = Joules / time of sat
    tSat = SatPulseNumberTotal*Params.pulseDur;
    Psat = Jsat / tSat;
    
    % Reorganize the Power equation to solve for B1 
    % B1 = sqrt(Psat/(epsilon*w0^2))
    % OLD WAY Params.satRMS = sqrt(Psat/(empFact*w0^2));
end

%% From here, calculate the flip angle based on the pulse shape and power
Psat = Psat/(empFact*w0^2); % rescale with empirical factor.

% Compute a temporary pulse to get the B1 integral from shape and duration
% Following Soustelle et al 2022, integral = p2, power = p2*B1peak^2
% Solve for peak, then scale normalized B1 and get integral for flipangle

if ~isfield(Params,'PulseOpt')
    switch Params.PulseOpt                
        % Special cases
        case 'gausshann'
            Params.PulseOpt.bw = 0.0002/Params.pulseDur; % override default Hann pulse shape.
        otherwise
            Params.PulseOpt = [];
    end
end

t = 0:Params.pulseDur/100:Params.pulseDur;
tempPulse = GetPulse(100, Params.delta, Params.pulseDur, ...
Params.SatPulseShape, Params.PulseOpt);

rf = tempPulse.('b1')(t);
p2 = trapz(t,rf.^2)/ Params.pulseDur; % See Soustelle et al 2022 for def.
B1peak = sqrt( Psat/p2);

if B1peak > B1peak_limit % conform to hardware constraint
    B1peak = B1peak_limit;
end

Params.satFlipAngle = trapz( t, rf*B1peak)*gam * 360;
Params.satB1peak = B1peak*1e6; % convert to microTesla


% figure; plot(t, rf*B1peak*1e6)

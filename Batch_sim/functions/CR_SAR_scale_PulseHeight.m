function Params = CR_SAR_scale_PulseHeight(Params)
% Uses SAR restrictions to give an estimate of the maximum pulse height 
% that can be used for the saturation pulses
% necessary parameters:
% numSat = number of sat pulses
% satRMS = root mean square of sat pulses (in Tesla)
% tp_sat = time of sat pulse (in seconds)
% numExc = number of excitation pulses
% flip = flip angle of excitation pulses (in degrees) 
% TR = time (in seconds)

if ~isfield(Params,'B0') % if not defined, assume 3T
    Params.B0 = 3; % main field strength (in Tesla)
end

if ~isfield(Params, 'boosted')
    error ('Please specify if boosted sat scheme is used (enter Params.boosted = 0 or 1)')
end

if Params.B0 == 7
    if strcmp(Params.TransmitCoil, 'STX')
        SAR_limit = 2.35; % Empirical value to match what I get at scanner
    end
else
    SAR_limit = 3; %(W/kg)
end
    
%% Multiply all variables except for the B1 field
% Using equation from Ibrahim, Tamer S. 2004. 
% Specify variables
r = 0.076; % meters for human head
l = 0.2; % meters for human head
w = 42.58e6 *Params.B0;

% There is likely time involved in the SAR monitoring for heat diffusion. A
% simplistic approach used here: two conductivity values to separate the
% protocols based on longer vs shorter TRs

if Params.B0 == 3
    % S/m from McCann et al 2019; 
    %conduc => % CSF = 1.71, GM is 0.466, WM is 0.21;
    if Params.boosted % modify for different definition of numSatPulse
        % conduc = 0.957; % empirical value TR = 3;
        % conduc = 0.81; % empirical value TR = 1.14
        conduc = 0.079*Params.TR + 0.72; % line fit to the above
    else
        conduc = 1.032; % empirical value TR = 0.12
    end
elseif Params.B0 == 7
    % From Van Lier et al 2014, the average conductivity at 7T was 
    % 13 percent higher than at 3T, when looking at GM,WM and CSF
    %conduc = 1.71*1.13;
    conduc = 1.48; % empirically solved based on what I get at scanner.
elseif Params.B0 == 1.5
    % From Van Lier et al 2014, the average conductivity at 1.5T was 
    % 10 percent lower than at 3T, when looking at GM,WM and CSF
    conduc = 1.71*0.9;
end

% calculate the power deposited
P_c = pi*r^4 * l * w^2 *conduc;

%% Power = J/s. Multiply by pulse time to find J of work done
% For Excitation pulse
% numberSatPulses *(powerCoefficient*B1field) * time pulse
excB1 = (Params.flipAngle*pi/180) / (42.58e6 * Params.WExcDur);
J_exc = Params.numExcitation*(P_c * excB1^2)* Params.WExcDur;

%% Combine the work, divide by TR
% Head SAR restriction is ~ 3 W/kg for head


if Params.boosted % modify for different definition of numSatPulse
    SatPulseNumberTotal = Params.numSatPulse * Params.satTrainPerBoost;
    Params.satRMS = abs(sqrt((SAR_limit*Params.TR - J_exc)/ (SatPulseNumberTotal*P_c*Params.pulseDur)));
else
    Params.satRMS = abs(sqrt((SAR_limit*Params.TR - J_exc)/(Params.numSatPulse*P_c*Params.pulseDur)));
end








if ~isfield(Params,'B0') % if not defined, assume 3T
    Params.B0 = 3; % main field strength (in Tesla)
end

if ~isfield(Params, 'boosted')
    error ('Please specify if boosted sat scheme is used (enter Params.boosted = 0 or 1)')
end


SAR_limit = 3; %(W/kg)
    
%% Multiply all variables except for the B1 field
% Using equation from Ibrahim, Tamer S. 2004. 
% Specify variables
w = 42.58e6 *Params.B0;

% There is likely time involved in the SAR monitoring for heat diffusion. A
% simplistic approach used here: two conductivity values to separate the
% protocols based on longer vs shorter TRs

if Params.B0 == 3
    r = 0.071; % meters for human head
    l = 0.57; % meters for human head

    % S/m from McCann et al 2019; 
    %conduc => % CSF = 1.71, GM is 0.466, WM is 0.21;
    if Params.boosted % modify for different definition of numSatPulse
        % conduc = 0.83; % empirical value TR = 3;
        % conduc = 0.805; % empirical value TR = 1.14
        conduc = 0.0134*Params.TR + 0.7897; % line fit to the above
    else
        conduc = 0.805; % empirical value TR = 0.100
        %conduc = 0.8; % empirical value TR = 0.12
    end
elseif Params.B0 == 7

    r = 0.066; % meters for human head
    l = 0.22; % meters for human head

    % From Van Lier et al 2014, the average conductivity at 7T was 
    % 13 percent higher than at 3T, when looking at GM,WM and CSF
    %conduc = 1.71*1.13;
    conduc = 0.8; % empirically solved based on what I get at scanner.
elseif Params.B0 == 1.5
    r = 0.1; % meters for human head
    l = 0.256; % meters for human head

    % From Van Lier et al 2014, the average conductivity at 1.5T was 
    % 10 percent lower than at 3T, when looking at GM,WM and CSF
    conduc = 1.71*0.9;
end

% calculate the power deposited
P_c = pi*r^4 * l * w^2 *conduc;

%% Power = J/s. Multiply by pulse time to find J of work done
% For Excitation pulse
% numberSatPulses *(powerCoefficient*B1field) * time pulse
excB1 = (Params.flipAngle) / (360*42.58e6 * Params.WExcDur);
J_exc = Params.numExcitation*(P_c * excB1^2)* Params.WExcDur;

%% Combine the work, divide by TR
% Head SAR restriction is ~ 3 W/kg for head
% Edge case that we use all the SAR up for excitation, check for negative
checkNeg = (SAR_limit*Params.TR - J_exc);

if checkNeg < 0
    Params.satRMS = -1;
else
    if Params.boosted % modify for different definition of numSatPulse
        SatPulseNumberTotal = Params.numSatPulse * Params.satTrainPerBoost;
        Params.satRMS = sqrt((SAR_limit*Params.TR - J_exc)/ (SatPulseNumberTotal*P_c*Params.pulseDur));
    else
        Params.satRMS = sqrt((SAR_limit*Params.TR - J_exc)/(Params.numSatPulse*P_c*Params.pulseDur));
    end
end


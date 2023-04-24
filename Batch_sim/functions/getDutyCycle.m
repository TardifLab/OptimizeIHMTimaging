function [DCsat, DCtr] = getDutyCycle(TR, TR_MT, numSatPulses, numExcitationPulses,...
        satPulseDur, excitationPulseDur, satPulseGap, Nburst)

% DC = duty cycle, is reported differently, so we will calculate the duty
% cycle of just the saturation period, and of the 

% using equation 2 from: A strategy to reduce the sensitivity of inhomogeneous
% magnetization transfer ( ihMT ) imaging to radiofrequency transmit field 
% variations at 3T. Soustelle et al., 2021. MRM.

% if only a single train of sat pulses are played out, TR_MT and Nburst
% might be previously undefined. 
if TR_MT == 0
    TR_MT = (satPulseDur+satPulseGap).*numSatPulses - satPulseGap; % remove one gap
end

if Nburst == 0
    Nburst = 1;
end

% in multitrain sat schemes, we drop the delay on the last one. Calculate
% for removal.
TD_MT = TR_MT - (satPulseDur+satPulseGap).*numSatPulses;
DCsat = (Nburst.*satPulseDur.*numSatPulses) ./ (Nburst.*TR_MT - TD_MT);

% Calculate time spent with RF on
tRF = excitationPulseDur.*numExcitationPulses + Nburst.*satPulseDur.*numSatPulses;
DCtr = tRF./TR;
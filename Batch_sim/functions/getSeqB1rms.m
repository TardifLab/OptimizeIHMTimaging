function [B1rmsSatOut, B1rmsTR] = getSeqB1rms(TR, TR_MT, numSatPulses, numExcitationPulses,...
        satPulseDur, excitationPulseDur, satPulseGap, Nburst, B1rmsSat, B1rmsExc, boosted)

% B1rms is reported differently, so we will calculate the b1rms of just the 
% saturation period, and of the entire TR 

% using equation 2 from: A strategy to reduce the sensitivity of inhomogeneous
% magnetization transfer ( ihMT ) imaging to radiofrequency transmit field 
% variations at 3T. Soustelle et al., 2021. MRM.


if TR_MT == 0
    TR_MT = (satPulseDur+satPulseGap).*numSatPulses;
end

if Nburst == 0
    Nburst = 1;
end

% in multitrain sat schemes, we drop the delay on the last one. Calculate
% for removal.
if boosted
    TD_MT = TR_MT - (satPulseDur+satPulseGap).*numSatPulses - satPulseGap;
else
    TD_MT = 0;
end

tSat = (Nburst.*TR_MT - TD_MT);

B1rmsSatOut = B1rmsSat.* sqrt( (Nburst.*satPulseDur.*numSatPulses) ./ tSat);

% Calculate B1rms of excitation period
tExc = (TR - tSat);
B1rms2 = B1rmsExc.* sqrt( (numExcitationPulses.*excitationPulseDur) ./ tExc);

% Time averaged B1rms-> sqrt of the squared b1 values multiplied by their
% percent contributions

B1rmsTR = sqrt( B1rmsSatOut.^2 .*tSat./TR + B1rms2.^2 .*tExc./TR);












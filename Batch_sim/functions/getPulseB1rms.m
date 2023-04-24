function B1rms = getPulseB1rms(satFlipAngle, pulseDur, SatPulseShape)

% This has been updated to take in the satFlipAngle, and then compute the
% B1rms using qMRlab pulse functions. To do so requires:
% 1. flip angle of the saturation pulse
% 2. duration of the saturation pulse (in seconds)
% 3. shape of the saturation pulse -> see qMRlab's 'GetPulse.m'

% The exported B1rms is in microTesla.


% Set up custom adjustments to pulses. Use with qMRlabs 'ViewPulse.m'
if strcmp(SatPulseShape, 'gausshann')
    PulseOpt.bw = 0.0002/pulseDur;
else
    PulseOpt = [];
end

delta = 5000; % not important for this. 
gam = 42.576;

tSat = 0:pulseDur/200:pulseDur;
satPulse = GetPulse(satFlipAngle, delta, pulseDur, SatPulseShape, PulseOpt);


B1_time = satPulse.omega(tSat)/(2*pi*gam);
P3 = (1/pulseDur) * trapz(tSat,B1_time.^2);
B1rms = sqrt(P3);




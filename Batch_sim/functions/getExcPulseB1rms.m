function B1rms = getExcPulseB1rms( flipAngle, pulseDur, PulseShape)

% I have only computed for one pulse shape, other can be done by going
% through the code commented out at the bottom.
% Currently assumes square (hard) pulse with 0.1ms duration

% using equation 3 from: A strategy to reduce the sensitivity of inhomogeneous magnetization transfer ( ihMT ) imaging to radiofrequency transmit field variations at 3 T
% requires qMRlab for getting values. 

% Pulse dur must be in seconds!


if ~exist('PulseShape','var')
    PulseShape = 'hard';
end

delta = 0; % not important for this. 
gam = 42.576; % MHz/T

tSat = 0:pulseDur/200:pulseDur;
excPulse = GetPulse(flipAngle, delta, pulseDur, PulseShape, []);


B1_time = excPulse.omega(tSat)/(2*pi*gam);
P3 = (1/pulseDur) * trapz(tSat,B1_time.^2);
B1rms = sqrt(P3);


%figure; ViewPulse(excPulse, 'omega')









function [Params, Params2, Params3] = getParamsForFitting()

% Used to clean up main function, have all the parameters set here. 
Params.B0 = 3;
Params.MTC = 1; % Magnetization Transfer Contrast
Params.TissueType = 'GM';
Params = DefaultCortexTissueParams(Params);
Params = CalcImagingParams(Params);


Params.b1 = 11.4; % microTesla
Params.numSatPulse = 6;
Params.pulseDur = 0.768/1000; %duration of 1 MT pulse in seconds
Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train
Params.TR = 120/1000; % total repetition time = MT pulse train and readout.
Params.WExcDur = 0.1/1000; % duration of water pulse
Params.numExcitation = 8; % number of readout lines/TR
Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
Params.delta = 8000;
Params.flipAngle = 5; % excitation flip angle water.
Params.echoSpacing = 7.66/1000;
Params.TD = Params.TR - (Params.numSatPulse *(Params.pulseDur+Params.pulseGapDur)) - (Params.numExcitation*Params.echoSpacing);
Params.SatPulseShape = 'gausshann';
Params.PulseOpt.bw = 0.3./Params.pulseDur; % override default Hann pulse shape.
Params.PerfectSpoiling = 1; % Speed up sims

Params.DummyEcho = 2;
Params.numExcitation = Params.numExcitation + Params.DummyEcho; % number of readout lines/TR WITH dummy

Params.boosted = 0; % use gaps in RF sat train -> NOTE this modifies the definition of numSatPul
Params.satTrainPerBoost = 1; % total number of pulses per TR = numSatPul *SatTrainPerBoost
Params.TR_MT = 0; % repetition time of satpulse train in seconds
Params = CalcVariableImagingParams(Params);


Params.NumLines = 216;
Params.NumPartitions = 192; 
Params.Slices = 176;
Params.Grappa = 1;
Params.ReferenceLines = 32;
Params.AccelerationFactor = 2;
Params.Segments = []; 
Params.TurboFactor = Params.numExcitation- Params.DummyEcho;
Params.ellipMask = 1;

%% Other 2 sequences:

Params2 = Params;
Params3 = Params;

Params2.b1 = 13.3; % microTesla
Params2.numSatPulse = 10;
Params2.pulseDur = 0.768/1000; %duration of 1 MT pulse in seconds
Params2.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train
Params2.TR = 3000/1000; % total repetition time = MT pulse train and readout.
Params2.WExcDur = 0.1/1000; % duration of water pulse
Params2.numExcitation = 200; % number of readout lines/TR
Params2.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
Params2.delta = 8000;
Params2.flipAngle = 11; % excitation flip angle water.
Params2.echoSpacing = 7.66/1000;
Params2.SatPulseShape = 'gausshann';
Params2.PulseOpt.bw = 0.3./Params2.pulseDur; % override default Hann pulse shape.

Params2.DummyEcho = 2;
Params2.numExcitation = Params2.numExcitation + Params2.DummyEcho; % number of readout lines/TR WITH dummy

Params2.boosted = 1; % use gaps in RF sat train -> NOTE this modifies the definition of numSatPul
Params2.satTrainPerBoost = 10; % total number of pulses per TR = numSatPul *SatTrainPerBoost
Params2.TR_MT = 90/1000; % repetition time of satpulse train in seconds
Params2.TD_MT = Params2.TR_MT - Params2.numSatPulse* (Params2.pulseDur + Params2.pulseGapDur) ;
Params2 = CalcVariableImagingParams(Params2);

Params2.TurboFactor = Params2.numExcitation- Params2.DummyEcho;

%%
Params3.b1 = 11.6; % microTesla
Params3.numSatPulse = 6;
Params3.pulseDur = 0.768/1000; %duration of 1 MT pulse in seconds
Params3.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train
Params3.TR = 1140/1000; % total repetition time = MT pulse train and readout.
Params3.WExcDur = 0.1/1000; % duration of water pulse
Params3.numExcitation = 80; % number of readout lines/TR
Params3.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
Params3.delta = 8000;
Params3.flipAngle = 7; % excitation flip angle water.
Params3.echoSpacing = 7.66/1000;
Params3.SatPulseShape = 'gausshann';
Params3.PulseOpt.bw = 0.3./Params3.pulseDur; % override default Hann pulse shape.

Params3.DummyEcho = 2;
Params3.numExcitation = Params3.numExcitation + Params3.DummyEcho; % number of readout lines/TR WITH dummy

Params3.boosted = 1; % use gaps in RF sat train -> NOTE this modifies the definition of numSatPul
Params3.satTrainPerBoost = 9; % total number of pulses per TR = numSatPul *SatTrainPerBoost
Params3.TR_MT = 60/1000; % repetition time of satpulse train in seconds
Params3 = CalcVariableImagingParams(Params3);

Params3.TD_MT = Params3.TR_MT - Params3.numSatPulse* (Params3.pulseDur + Params3.pulseGapDur) ;
Params3.TurboFactor = Params3.numExcitation- Params3.DummyEcho;




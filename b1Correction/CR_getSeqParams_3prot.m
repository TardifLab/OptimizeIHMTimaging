function [Params, outputSamplingTable] = CR_getSeqParams_3prot(turbofactor)

%filled with values from simulation results:
Params.B0 = 3;
Params.MTC = 1; % Magnetization Transfer Contrast
Params.TissueType = 'GM';
Params = DefaultCortexTissueParams(Params);
Params = CalcImagingParams(Params);

Params.WExcDur = 0.1/1000; % duration of water pulse
Params.echospacing = 7.66/1000;
Params.PerfectSpoiling = 1;

%% MT params
if  turbofactor == 8
    Params.delta = 8000;
    Params.flipAngle = 5; % excitation flip angle water.
    Params.TR = 0.12; % total repetition time = MT pulse train and readout.
    Params.numSatPulse = 6;
    Params.TurboFactor = 8;
    Params.pulseDur = 0.768/1000; %duration of 1 MT pulse in seconds
    Params.b1 = 11.4; % microTesla
    Params.satTrainPerBoost = 1;
    Params.TR_MT = 0; 
    Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
    Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train % C.R. new, shift from 1ms to 0.5
    Params.DummyEcho = 2;
    Params.numExcitation = Params.TurboFactor + Params.DummyEcho; % number of readout lines/TR
    Params.boosted = 0;
    Params.satTrainPerBoost = 1; 
    Params.TR_MT = 0; 
    
elseif turbofactor == 80
    Params.delta = 8000;
    Params.flipAngle = 7; % excitation flip angle water.
    Params.TR = 1.14; % total repetition time = MT pulse train and readout.
    Params.numSatPulse = 6;
    Params.TurboFactor = 80;
    Params.pulseDur = 0.768/1000; %duration of 1 MT pulse in seconds
    Params.b1 = 11.6; % microTesla
    Params.satTrainPerBoost = 9;
    Params.TR_MT = 0.06;    
    Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
    Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train % C.R. new, shift from 1ms to 0.5
    Params.DummyEcho = 2;
    Params.numExcitation = Params.TurboFactor + Params.DummyEcho; % number of readout lines/TR
    Params.boosted = 1;

elseif turbofactor == 200
    Params.delta = 8000;
    Params.flipAngle = 11; % excitation flip angle water.
    Params.TR = 3.0; % total repetition time = MT pulse train and readout.
    Params.numSatPulse = 10;
    Params.TurboFactor = 200;
    Params.pulseDur = 0.768/1000; %duration of 1 MT pulse in seconds
    Params.b1 = 13.3; % microTesla
    Params.satTrainPerBoost = 10;
    Params.TR_MT = 90/1000;
    Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
    Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train % C.R. new, shift from 1ms to 0.5
    Params.DummyEcho = 2;
    Params.numExcitation = Params.TurboFactor + Params.DummyEcho; % number of readout lines/TR
    Params.boosted = 1;

else
    error('turbofactor input was not one of the options')

end 

% MT parameters that will be consistent:
Params.SatPulseShape = 'gausshann';
Params.PulseOpt.bw = 0.3./Params.pulseDur; % override default Hann pulse shape.
Params.TD_MT =  Params.TR_MT - Params.numSatPulse* (Params.pulseDur + Params.pulseGapDur) ;   
Params = CalcVariableImagingParams(Params);


% Other image parameters
Params.NumLines = 216;
Params.NumPartitions = 192; 
Params.Slices = 176;
Params.Grappa = 1;
Params.ReferenceLines = 32;
Params.AccelerationFactor = 2;
Params.Segments = []; 
Params.TurboFactor = Params.numExcitation- Params.DummyEcho;
Params.ellipMask = 1;

[outputSamplingTable, ~, Params.Segments] = Step1_calculateKspaceSampling_v3 (Params);




































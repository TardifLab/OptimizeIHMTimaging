function Params = CR_getSeqParams(turbofactor)

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
if  turbofactor == 200
    Params.delta = 8000;
    Params.flipAngle = 11; % excitation flip angle water.
    Params.TR = 2.9; % total repetition time = MT pulse train and readout.
    Params.numSatPulse = 8;
    Params.TurboFactor = 200;
    Params.pulseDur = 0.768/1000; %duration of 1 MT pulse in seconds
    Params.b1 = 13.9; % microTesla
    Params.satTrainPerBoost = 12;
    Params.TR_MT = 0.04; 
    Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
    Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train % C.R. new, shift from 1ms to 0.5
    Params.DummyEcho = 2;
    Params.numExcitation = Params.TurboFactor + Params.DummyEcho; % number of readout lines/TR
    Params.boosted = 1;

elseif turbofactor == 160
    Params.delta = 8000;
    Params.flipAngle = 11; % excitation flip angle water.
    Params.TR = 2.3; % total repetition time = MT pulse train and readout.
    Params.numSatPulse = 6;
    Params.TurboFactor = 160;
    Params.pulseDur = 0.768/1000; %duration of 1 MT pulse in seconds
    Params.b1 = 13.8; % microTesla
    Params.satTrainPerBoost = 12;
    Params.TR_MT = 0.04;    
    Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
    Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train % C.R. new, shift from 1ms to 0.5
    Params.DummyEcho = 2;
    Params.numExcitation = Params.TurboFactor + Params.DummyEcho; % number of readout lines/TR
    Params.boosted = 1;

elseif turbofactor == 120
    Params.delta = 8000;
    Params.flipAngle = 5; % excitation flip angle water.
    Params.TR = 1.75; % total repetition time = MT pulse train and readout.
    Params.numSatPulse = 6;
    Params.TurboFactor = 120;
    Params.pulseDur = 0.768/1000; %duration of 1 MT pulse in seconds
    Params.b1 = 14; % microTesla
    Params.satTrainPerBoost = 14;
    Params.TR_MT = 0.04;
    Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
    Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train % C.R. new, shift from 1ms to 0.5
    Params.DummyEcho = 2;
    Params.numExcitation = Params.TurboFactor + Params.DummyEcho; % number of readout lines/TR
    Params.boosted = 1;

elseif turbofactor == 90
    Params.delta = 8000;
    Params.flipAngle = 5; % excitation flip angle water.
    Params.TR = 1.25; % total repetition time = MT pulse train and readout.
    Params.numSatPulse = 6;
    Params.TurboFactor = 90;
    Params.pulseDur = 0.768/1000; %duration of 1 MT pulse in seconds
    Params.b1 = 13.7; % microTesla
    Params.satTrainPerBoost = 11;
    Params.TR_MT = 0.04;
    Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
    Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train % C.R. new, shift from 1ms to 0.5
    Params.DummyEcho = 2;
    Params.numExcitation = Params.TurboFactor + Params.DummyEcho; % number of readout lines/TR
    Params.boosted = 1;

elseif turbofactor == 80
    Params.delta = 8000;
    Params.flipAngle = 5; % excitation flip angle water.
    Params.TR = 1.25; % total repetition time = MT pulse train and readout.
    Params.numSatPulse = 6;
    Params.TurboFactor = 80;
    Params.pulseDur = 0.768/1000; %duration of 1 MT pulse in seconds
    Params.b1 = 14; % microTesla
    Params.satTrainPerBoost = 11;
    Params.TR_MT = 0.04;
    Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
    Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train % C.R. new, shift from 1ms to 0.5
    Params.DummyEcho = 2;
    Params.numExcitation = Params.TurboFactor + Params.DummyEcho; % number of readout lines/TR
    Params.boosted = 1;

elseif turbofactor == 48
    Params.delta = 8000;
    Params.flipAngle = 5; % excitation flip angle water.
    Params.TR = 0.75; % total repetition time = MT pulse train and readout.
    Params.numSatPulse = 6;
    Params.TurboFactor = 48;
    Params.pulseDur = 0.768/1000; %duration of 1 MT pulse in seconds
    Params.b1 = 14; % microTesla
    Params.satTrainPerBoost = 7;
    Params.TR_MT = 0.04;
    Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
    Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train % C.R. new, shift from 1ms to 0.5
    Params.DummyEcho = 2;
    Params.numExcitation = Params.TurboFactor + Params.DummyEcho; % number of readout lines/TR
    Params.boosted = 1;
    
elseif turbofactor == 10
    Params.delta = 8000;
    Params.flipAngle = 4; % excitation flip angle water.
    Params.TR = 0.140; % total repetition time = MT pulse train and readout.
    Params.numSatPulse = 6;
    Params.TurboFactor = 10;
    Params.pulseDur = 0.768/1000; %duration of 1 MT pulse in seconds
    Params.b1 = 13.7; % microTesla
    Params.satTrainPerBoost = 0;
    Params.TR_MT = 0; 
    Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
    Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train % C.R. new, shift from 1ms to 0.5
    Params.DummyEcho = 2;
    Params.numExcitation = Params.TurboFactor + Params.DummyEcho; % number of readout lines/TR
    Params.boosted = 0;
else
    error('turbofactor input was not one of the options')

end 

% MT parameters that will be consistent:
Params.SatPulseShape = 'gausshann';
Params.PulseOpt.bw = 0.3./Params.pulseDur; % override default Hann pulse shape.
Params.TD_MT =  Params.TR_MT - Params.numSatPulse* (Params.pulseDur + Params.pulseGapDur) ;   


Params = CalcVariableImagingParams(Params);




%% Pull parameters using:
% tempS = simResults;
% tempP = convParam;
% 
% TF2 =  tempP(:,6) ~= 10; 
% tempS(TF2,:) = [];
% tempP(TF2,:) = [];
% 
% 
% [temp, top10Effidx] = sort(tempS(:,11), 'descend'); % sort by ihMTsat SNR efficiency
% Top10sorted_b = tempP ( top10Effidx,:);
% Top10sorted_b = array2table(Top10sorted_b, 'VariableNames',{'delta', 'flipAngle', ...
%     'TR', 'numSatPulse', 'pulseDur', 'numExc', 'satTrainPerBoost', 'TR_MT','b1'});
% 
% 
% snr_t = temp(1)































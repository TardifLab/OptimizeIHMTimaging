%% Goal of this script is to simulate many options for generating ihMT contrast
% to determine which is optimal for the cortex.

% Version 4 is a complete re-write to optimize with parpool on a cluster with many cores available

% Version 3 adds the ability to go for a boosted approach
%
% Version 2 adds dummy scans to mimic the TFL sequence. And columns for MTsat

%% After this, generate figures with plotSimResultsIHMT_figures_fineGridDummyv3.m

% location of 
baseDir = '/data_/tardiflab/chris/development/scripts/Deichman2000_backup/cortical_ihMT_sim/sim_wSpoil/RF_grad_diffusion_v4/';

addpath(genpath(baseDir)) % new sim code.
addpath(genpath('/data_/tardiflab/chris/development/scripts/Deichman2000_backup/cortical_ihMT_sim/sim_wSpoil/qMRLab-master')) %qMRlab
  

DATADIR = strcat(baseDir,'Batch_sim/3T_v1/outputsConventional/');
DATADIR_INT = strcat(DATADIR,'intermediate/');
mkdir(DATADIR)
mkdir(DATADIR_INT)

% Import k-space files 
ref_kspace_dir = strcat(baseDir,'kspaceWeighting/Atlas_reference/');
load( strcat( ref_kspace_dir,'GM_seg_MNI_152_image.mat')) % keep these here so you know what files you are setting directory to
load( strcat( ref_kspace_dir,'GM_seg_MNI_152_kspace.mat'))


cl = parcluster('local');
cl.NumWorkers = 20;
saveProfile(cl);

parpool(str2num(getenv('NSLOTS')));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF MODIFY FOR CUSTOM PATHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up Param Structure.
Params.B0 = 3;
Params.MTC = 1; % Magnetization Transfer Contrast
Params.TissueType = 'GM';
Params = DefaultCortexTissueParams(Params);
Params = CalcImagingParams(Params);


%% Custom fit cortex parameters, should already be set in DefaultCortexTissueParams()



%% Sequence Parameters
Params.B0 = 3; % in Tesla
Params.b1 = []; % microTesla
Params.numSatPulse = [];
Params.pulseDur = []; %duration of 1 MT pulse in seconds
Params.TR = []; % total repetition time = MT pulse train and readout.
Params.numExcitation = []; % number of readout lines/TR
Params.flipAngle = []; % excitation flip angle water.
Params.delta = [];
Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train % C.R. new, shift from 1ms to 0.5
Params.WExcDur = 0.1/1000; % duration of water pulse
Params.echoSpacing = 7.66/1000;
Params.ReferenceScan = 0;
Params.SatPulseShape = 'gausshann';
% Params.PulseOpt.bw = 0.3./Params.pulseDur; % override default Hann pulse shape.
Params.PerfectSpoiling = 1;

Params = CalcVariableImagingParams(Params);


%% Acquistion matrix 
Params.NumLines = 216;
Params.NumPartitions = 192; 
Params.Slices = 176; 
Params.Grappa = 1;
Params.ReferenceLines = 32;
Params.AccelerationFactor = 2;
Params.Segments = []; 
Params.TurboFactor = []; %Params.numExcitation- Params.DummyEcho;
Params.ellipMask = 1;


%% Simulation Variable

Offset_freq = 8000; 
flipAngle = [4,5,6, 7, 8, 9, 10, 11, 13];
pulseDur = [0.768, 1.024,1.28, 2.048]./1000; % in seconds
TR = [15, 30, 60,70, 80,90, 100, 110, 120, 140, 180, 250, 500, 1140]./1000; % in seconds
numSatPulse = 2:2:20;
numExc = [1, 2, 4, 6, 8, 10,12,14, 16,18, 20,22, 24, 28,32, 48, 64, 90];
Params.boosted = 0; % 
Params.satTrainPerBoost = 0;
satTrainPerBoost = 1; % Not applicable here
TR_MT = 1; % Not applicable here
Params.TR_MT = 0;
B1rms_limit = 14e-6; % in Tesla

c1 = length(Offset_freq);
c2 = length(flipAngle);
c3 = length(TR);
c4 = length(numSatPulse);
c5 = length(pulseDur);
c6 = length(numExc);
c7 = length(satTrainPerBoost);
c8 = length(TR_MT);


simLength = c1*c2*c3*c4*c5*c6*c7*c8;


%% First make parameter set variable that stores sequence options.
ParameterSet = zeros(simLength,8);
idx = 1;

% cutoff time
imgTime = 5*60; 
refVal = imgTime/((Params.NumLines/Params.AccelerationFactor+Params.ReferenceLines) * Params.NumPartitions *0.7) ; 


for a = 1:c1
    for b = 1:c2
        for c = 1:c3
            for d = 1:c4
                for e = 1:c5
                    for f = 1:c6
                        for g = 1:c7
                            for h = 1:c8

                                if numExc(f) == 1
                                   DummyEcho = 0;
                                elseif numExc(f) < 4
                                    DummyEcho = 1;
                                else 
                                    DummyEcho = 2;
                                end

                                chkTime = numSatPulse(d)*(pulseDur(e)+Params.pulseGapDur)+...
                                    (DummyEcho+numExc(f))*(Params.echoSpacing)- ...
                                    Params.pulseGapDur + Params.G_time_elapse_MT; % last term is dephasing gradient

                                if (TR(c)/numExc(f) < refVal) && (chkTime < TR(c)) % then parameters will fit within scan time
                                    ParameterSet(idx,:) = [Offset_freq(a),...
                                        flipAngle(b), TR(c), numSatPulse(d),...
                                        pulseDur(e), numExc(f), satTrainPerBoost(g),...
                                        TR_MT(h) ]; 
        
                                    idx = idx+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

ParameterSet(idx:end,:) = [];
simLength = idx-1;
save( strcat(DATADIR,'ParameterSet_3T_batch_conventional.mat'),'ParameterSet')



%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECTION SHOULD BE CORRECT %%%%%%%%%%%%%%%%%%%%%%%
% Note a hack was needed to get around parfor restrictions:
%     reference_Signal(qi,1:ParameterSet(qi,6)) = referenceSignal;

% Turned into this:
%     tV = zeros(1,maxExc);
%     tV(1:ParameterSet(qi,6)) = referenceSignal;
%     reference_Signal(qi,:) = tV;


%% Loop through ParameterSet to get magnetization vector over excitation train

maxExc = max(numExc);

single_Signal = zeros(simLength, maxExc);
single_B1_val = zeros(simLength, 1);
dual_Signal = zeros(simLength, maxExc);
dual_B1_val = zeros(simLength, 1);
reference_Signal = zeros(simLength, maxExc);


% Need to break into chunks to allow for saving.
loopVec = 1:simLength;
numChunk = 10;
st = round(linspace( min(loopVec),max(loopVec), numChunk +1));
ed = st(2:end); % remove start index
ed(1:end-1) = ed(1:end-1)-1; % increment to prevent overlap
st = st(1:end-1); % remove last index


tic
for i = 1:numChunk
    stVal = st(i);
    edVal = ed(i);

    parfor qi = stVal:edVal
    
        % make calls to BlochSimFlashSequence_v1
        [inputMag, B1_val] = CR_batch_simSequenceFunction(Params, ...
            'delta', ParameterSet(qi,1),...
            'flipAngle', ParameterSet(qi,2),...
            'TR', ParameterSet(qi,3),...
            'numSatPulse', ParameterSet(qi,4),...
            'pulseDur', ParameterSet(qi,5),...
            'numExcitation', ParameterSet(qi,6),...
            'satTrainPerBoost', ParameterSet(qi,7),...
            'TR_MT', ParameterSet(qi,8),...
            'freqPattern','single');
    
        tV = zeros(1,maxExc);
        tV(1:ParameterSet(qi,6)) = inputMag;
        single_Signal(qi,:) = tV;
        single_B1_val(qi) = B1_val;
    
        [inputMag, B1_val] = CR_batch_simSequenceFunction(Params, ...
            'delta', ParameterSet(qi,1),...
            'flipAngle', ParameterSet(qi,2),...
            'TR', ParameterSet(qi,3),...
            'numSatPulse', ParameterSet(qi,4),...
            'pulseDur', ParameterSet(qi,5),...
            'numExcitation', ParameterSet(qi,6),...
            'satTrainPerBoost', ParameterSet(qi,7),...
            'TR_MT', ParameterSet(qi,8),...
            'freqPattern','dualAlternate');
    
        tV = zeros(1,maxExc);
        tV(1:ParameterSet(qi,6)) = inputMag;
        dual_Signal(qi,:) = tV;
        dual_B1_val(qi) = B1_val;
    
        
        % If at least one of the signal values is >0, calculate reference signal for MTR
        if max(single_Signal(qi,:)) > 0 
            [referenceSignal, ~] = CR_batch_simSequenceFunction(Params, ...
                'delta', ParameterSet(qi,1),...
                'flipAngle', ParameterSet(qi,2),...
                'TR', ParameterSet(qi,3),...
                'numSatPulse', ParameterSet(qi,4),...
                'pulseDur', ParameterSet(qi,5),...
                'numExcitation', ParameterSet(qi,6),...
                'satTrainPerBoost', ParameterSet(qi,7),...
                'TR_MT', ParameterSet(qi,8),...
                'freqPattern','dualAlternate', ...
                'ReferenceScan', 1); % add to calculate MTR
        else 
            referenceSignal = 0;
        end
    
        tV = zeros(1,maxExc);
        tV(1:ParameterSet(qi,6)) = referenceSignal;
        reference_Signal(qi,:) = tV;
    
    end

    % In the event that things can go wrong, lets periodically save intermediate results
    save( strcat(DATADIR_INT,'nonBoost_reference_Signal_intermed_',num2str(i),'.mat'),'reference_Signal')
    save( strcat(DATADIR_INT,'nonBoost_dual_Signal_intermed_',num2str(i),'.mat'),'dual_Signal')
    save( strcat(DATADIR_INT,'nonBoost_dual_B1_val_intermed_',num2str(i),'.mat'),'dual_B1_val')
    save( strcat(DATADIR_INT,'nonBoost_single_Signal_intermed_',num2str(i),'.mat'),'single_Signal')
    save( strcat(DATADIR_INT,'nonBoost_single_B1_val_intermed_',num2str(i),'.mat'),'single_B1_val')

    disp( i/numChunk *100) % print percent done
    toc

end
toc % DONE SIMULATIONS
          
if max(dual_B1_val) == 0
    error('Didnt work')
end

save( strcat(DATADIR,'nonBoost_reference_Signal.mat'),'reference_Signal')
save( strcat(DATADIR,'nonBoost_dual_Signal.mat'),'dual_Signal')
save( strcat(DATADIR,'nonBoost_dual_B1_val.mat'),'dual_B1_val')
save( strcat(DATADIR,'nonBoost_single_Signal.mat'),'single_Signal')
save( strcat(DATADIR,'nonBoost_single_B1_val.mat'),'single_B1_val')






%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECTION SHOULD BE CORRECT %%%%%%%%%%%%%%%%%%%%%%%
%% Loop through the Magnetization vectors to extract the grey matter signal


single_GM_val = zeros(simLength,1);
dual_GM_val = zeros(simLength,1);
ref_GM_val  = zeros(simLength,1);

% speed up by sorting based on number excitation, then only computing the sampling table once per turbofactor
for i = 1:c6 % length(numExc);

    nE = numExc(i);

    % Calculate sampling table for the selected turbofactor
    Params.TurboFactor = nE; % - DummyEcho;
    [outputSamplingTable, ~, Params.Segments] = Step1_calculateKspaceSampling_v3 (Params);

    % Find entries with the selected turbofactor
    q = find( (ParameterSet(:,6) == nE ) );
    tempParams = ParameterSet(q,:);

    tempSingleSig = zeros( length(tempParams), 1); 
    tempDualSig   = zeros( length(tempParams), 1); 
    tempRefSig    = zeros( length(tempParams), 1); 

    disp(strcat('For Turbofactor = ', num2str(nE),'...' ))

    % Parallel loop over values
    parfor qi = 1:length(tempParams)

        % Extract signal values over the excitation train
        inputS = single_Signal(q(qi),1:nE);
        inputD = dual_Signal(q(qi),1:nE);
        inputR = reference_Signal(q(qi),1:nE);

        % Params details passed along are just related to image matrix size, so it should be the same between
        % all of the protocols
        tempSingleSig(qi) = CR_generate_BSF_scaling_v1( inputS, Params, outputSamplingTable, gm_m, fft_gm_m);
        tempDualSig(qi)   = CR_generate_BSF_scaling_v1( inputD, Params, outputSamplingTable, gm_m, fft_gm_m);
        tempRefSig(qi)    = CR_generate_BSF_scaling_v1( inputR, Params, outputSamplingTable, gm_m, fft_gm_m);

    end

    % Return the subsampled signal set to the larger matrix containing all turbofactors
    single_GM_val(q) = tempSingleSig;
    dual_GM_val(q)   = tempDualSig;
    ref_GM_val(q)    = tempRefSig;

    disp(strcat(num2str(i/c6 *100), '% done...' ))
end

save( strcat(DATADIR,'nonBoost_single_GM_val.mat'), 'single_GM_val')
save( strcat(DATADIR,'nonBoost_dual_GM_val.mat'),   'dual_GM_val')
save( strcat(DATADIR,'nonBoost_ref_GM_val.mat'),    'ref_GM_val')




%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECTION SHOULD BE CORRECT %%%%%%%%%%%%%%%%%%%%%%%
%% With the signal sorted out, then calculate MTsat, and ihMTsat
single_Sat_val = zeros(simLength,1);
dual_Sat_val = zeros(simLength,1);

% Pull out a few things for faster running
T1 = (1/Params.Raobs)*1000;
M0 = Params.M0a;
echoSpacing = Params.echoSpacing*1000;

parfor qi = 1: simLength

    if ParameterSet(qi,6) == 1
        DummyEcho = 0;
    elseif ParameterSet(qi,6) < 4
        DummyEcho = 1;
    else 
        DummyEcho = 2;
    end
    
    single_Sat_val(qi) =calcMTsatThruLookupTablewithDummyV3( single_GM_val(qi), ...
        [], T1, 1, M0, echoSpacing, ParameterSet(qi,6) + DummyEcho,...
        ParameterSet(qi,3)*1000, ParameterSet(qi,2), DummyEcho);  

    dual_Sat_val(qi) =calcMTsatThruLookupTablewithDummyV3( dual_GM_val(qi), ...
        [], T1, 1, M0, echoSpacing, ParameterSet(qi,6) + DummyEcho, ...
        ParameterSet(qi,3)*1000, ParameterSet(qi,2), DummyEcho);  

end

save( strcat(DATADIR,'nonBoost_single_Sat_val.mat'), 'single_Sat_val')
save( strcat(DATADIR,'nonBoost_dual_Sat_val.mat'),   'dual_Sat_val')






%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECTION NEEDS FILLED IN %%%%%%%%%%%%%%%%%%%%%%%
% Finish by calculating MTR, ihMTR and ihMTsat values. Combine into one matrix and save.


ihMT_sat=  dual_Sat_val - single_Sat_val ; 
sMTR  = (ref_GM_val - single_GM_val)./ref_GM_val;
dMTR  = (ref_GM_val - dual_GM_val)./ref_GM_val;
ihMTR = (single_GM_val - dual_GM_val)./ ref_GM_val; 



simResults = [single_GM_val, dual_GM_val, ref_GM_val, sMTR, dMTR, ihMTR,...
     single_Sat_val,dual_Sat_val, ihMT_sat];


save( strcat(DATADIR,'nonBoost_simResults_calcMetrics.mat'), 'simResults')



% 
% ihMT_proxy = (single_GM_val - dual_GM_val);
% max(ihMT_sat)





















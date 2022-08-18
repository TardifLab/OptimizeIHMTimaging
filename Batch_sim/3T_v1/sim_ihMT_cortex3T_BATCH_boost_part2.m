%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Goal of this script is to simulate many options for generating ihMT contrast
% to determine which is optimal for the cortex.


%% After this, generate figures with plotSimResultsIHMT_figures_fineGridDummyv3.m

baseDir = '/path/to/OptimizeIHMTimaging/';
addpath(genpath(baseDir)) % new sim code.


DATADIR = strcat(baseDir,'Batch_sim/3T_v1/outputsBoost/');

% Import k-space files 
ref_kspace_dir = strcat(baseDir,'kspaceWeighting/Atlas_reference/');
load( strcat( ref_kspace_dir,'GM_seg_MNI_152_image.mat')) % keep these here so you know what files you are setting directory to
load( strcat( ref_kspace_dir,'GM_seg_MNI_152_kspace.mat'))


%% This is heavier on memory than the previous script. 
% Depending on your RAM, you might need to lower the number of cores
% compared to PART1.
cl = parcluster('local');
cl.NumWorkers = 10;
saveProfile(cl);


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
Params.freqPattern = 'single'; % options: 'single', 'dualAlternate', 'dualContinuous'
Params.pulseGapDur = 0.3/1000; %ms gap between MT pulses in train
Params.WExcDur = 0.1/1000; % duration of water pulse
Params.echoSpacing = 7.66/1000;
Params.ReferenceScan = 0;
Params.SatPulseShape = 'gausshann';
% Params.PulseOpt.bw = 0.3./Params.pulseDur; % override default Hann pulse shape.

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

% For Munsch style number excitations sim
Offset_freq = 8000; %5000:1000:13000;
flipAngle = [5, 7, 8 9,10, 11] ; 
TR = [200, 350, 500, 750, 1000, 1250, 1500, 1750, 2000, 2300, 2600, 2900, 3200, 3600, 4000]./1000; % in seconds
pulseDur = [0.768, 1.024]./1000; % in seconds
numExc = [16, 24, 32, 48, 64, 80, 90, 100, 110, 120, 140, 160, 180, 200];
Params.boosted = 1; % use gaps in RF sat train -> NOTE this modifies the definition of numSatPul
numSatPulse = 4:2:20; % this is PER RF train (not necessarily in the whole TR )
satTrainPerBoost = 2:14; % total number of pulses per TR = numSatPul *SatTrainPerBoost
TR_MT = [40, 60,70, 80, 90, 100, 110, 120, 150]./1000; % repetition time of satpulse train in seconds

B1rms_limit = 14e-6; % in Tesla

c1 = length(Offset_freq);
c2 = length(flipAngle);
c3 = length(TR);
c4 = length(numSatPulse);
c5 = length(pulseDur);
c6 = length(numExc);
c7 = length(satTrainPerBoost);
c8 = length(TR_MT);

%% Load Previous Variables

% Parameter set.
DATADIR_INT = strcat(DATADIR,'intermediate/');
DATADIR1 = strcat(baseDir,'Batch_sim/3T_v1/outputsBoost/');
ParameterSet = load( strcat(DATADIR1,'ParameterSet_3T_batch_boost.mat'));
ParameterSet = ParameterSet.ParameterSet;

simLength = length(ParameterSet);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% You may need to adjust this. This assumed the PART1 script was run in 3 chunks


% Part1_1
 rs1 = load( strcat(DATADIR_INT,'Boost_reference_Signal_intermed_16.mat'),'reference_Signal'); 
 ds1 = load( strcat(DATADIR_INT,'Boost_dual_Signal_intermed_16.mat'),'dual_Signal');
 db1 = load( strcat(DATADIR_INT,'Boost_dual_B1_val_intermed_16.mat'),'dual_B1_val');
 ss1 = load( strcat(DATADIR_INT,'Boost_single_Signal_intermed_16.mat'),'single_Signal');
 sb1 = load( strcat(DATADIR_INT,'Boost_single_B1_val_intermed_16.mat'),'single_B1_val');

 rs1 = rs1.reference_Signal;
 ds1 = ds1.dual_Signal;
 db1 = db1.dual_B1_val;
 ss1 = ss1.single_Signal;
 sb1 = sb1.single_B1_val;
  
% Part1_2
 rs2 = load( strcat(DATADIR_INT,'Boost_reference_Signal_intermed_33.mat'),'reference_Signal'); 
 ds2 = load( strcat(DATADIR_INT,'Boost_dual_Signal_intermed_33.mat'),'dual_Signal');
 db2 = load( strcat(DATADIR_INT,'Boost_dual_B1_val_intermed_33.mat'),'dual_B1_val');
 ss2 = load( strcat(DATADIR_INT,'Boost_single_Signal_intermed_33.mat'),'single_Signal');
 sb2 = load( strcat(DATADIR_INT,'Boost_single_B1_val_intermed_33.mat'),'single_B1_val');

 rs2 = rs2.reference_Signal;
 ds2 = ds2.dual_Signal;
 db2 = db2.dual_B1_val;
 ss2 = ss2.single_Signal;
 sb2 = sb2.single_B1_val;

% Part1_3
 rs3 = load( strcat(DATADIR_INT,'Boost_reference_Signal_intermed_50.mat'),'reference_Signal'); 
 ds3 = load( strcat(DATADIR_INT,'Boost_dual_Signal_intermed_50.mat'),'dual_Signal');
 db3 = load( strcat(DATADIR_INT,'Boost_dual_B1_val_intermed_50.mat'),'dual_B1_val');
 ss3 = load( strcat(DATADIR_INT,'Boost_single_Signal_intermed_50.mat'),'single_Signal');
 sb3 = load( strcat(DATADIR_INT,'Boost_single_B1_val_intermed_50.mat'),'single_B1_val');

 rs3 = rs3.reference_Signal;
 ds3 = ds3.dual_Signal;
 db3 = db3.dual_B1_val;
 ss3 = ss3.single_Signal;
 sb3 = sb3.single_B1_val;

 
 %% Concatenate them:

loopVec = 1:simLength;
numChunk = 50;
st = round(linspace( min(loopVec),max(loopVec), numChunk +1));
ed = st(2:end); % remove start index
ed(1:end-1) = ed(1:end-1)-1; % increment to prevent overlap
st = st(1:end-1); % remove last index
 
reference_Signal = vertcat(rs1(st(1):ed(16),:),rs2(st(17):ed(33),:), rs3(st(34):ed(50),:));
dual_Signal = vertcat(ds1(st(1):ed(16),:),ds2(st(17):ed(33),:), ds3(st(34):ed(50),:));
dual_B1_val = vertcat(db1(st(1):ed(16),:),db2(st(17):ed(33),:), db3(st(34):ed(50),:));
single_Signal = vertcat(ss1(st(1):ed(16),:),ss2(st(17):ed(33),:), ss3(st(34):ed(50),:));
single_B1_val = vertcat(sb1(st(1):ed(16),:),sb2(st(17):ed(33),:), sb3(st(34):ed(50),:));
 

clear rs1 rs2 rs3 ds1 ds2 ds3 db1 db2 db3 ss1 ss2 ss3 sb1 sb2 sb3

save( strcat(DATADIR,'Boost_reference_Signal.mat'),'reference_Signal')
save( strcat(DATADIR,'Boost_dual_Signal.mat'),'dual_Signal')
save( strcat(DATADIR,'Boost_dual_B1_val.mat'),'dual_B1_val')
save( strcat(DATADIR,'Boost_single_Signal.mat'),'single_Signal')
save( strcat(DATADIR,'Boost_single_B1_val.mat'),'single_B1_val')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% 
% 
% load( strcat(DATADIR,'Boost_reference_Signal.mat'))
% load( strcat(DATADIR,'Boost_dual_Signal.mat'))
% load( strcat(DATADIR,'Boost_dual_B1_val.mat'))
% load( strcat(DATADIR,'Boost_single_Signal.mat'))
% load( strcat(DATADIR,'Boost_single_B1_val.mat'))

disp('Intermediate files loaded and concatenated')

%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%
%% Loop through the Magnetization vectors to extract the grey matter signal

parpool(str2num(getenv('NSLOTS')));

maxExc = max(numExc);

single_GM_val = zeros(simLength,1);
dual_GM_val = zeros(simLength,1);
ref_GM_val  = zeros(simLength,1);

% speed up by sorting based on number excitation, then only computing the sampling table once per turbofactor
disp('Running signal scaling function...')
for i = 1:c6 % length(numExc);

    nE = numExc(i);
   
    % Calculate sampling table for the selected turbofactor
    Params.TurboFactor = nE; 
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

save( strcat(DATADIR,'Boost_single_GM_val.mat'), 'single_GM_val')
save( strcat(DATADIR,'Boost_dual_GM_val.mat'),   'dual_GM_val')
save( strcat(DATADIR,'Boost_ref_GM_val.mat'),    'ref_GM_val')






%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECTION SHOULD BE CORRECT %%%%%%%%%%%%%%%%%%%%%%%
%% With the signal sorted out, then calculate MTsat, and ihMTsat
single_Sat_val = zeros(simLength,1);
dual_Sat_val = zeros(simLength,1);

% Pull out a few things for faster running
T1 = (1/Params.Raobs)*1000;
M0 = Params.M0a;
echoSpacing = Params.echoSpacing*1000;

disp('Now calculating MTsat' )
parfor qi = 1: simLength

    if ParameterSet(qi,6) == 1
        DummyEcho = 0;
    elseif ParameterSet(qi,6) < 4
        DummyEcho = 1;
    else 
        DummyEcho = 2;
    end
    
    single_Sat_val(qi) =calcMTsatThruLookupTablewithDummyV3( single_GM_val(qi), ...
        1, T1, 1, M0, echoSpacing, ParameterSet(qi,6) + DummyEcho,...
        ParameterSet(qi,3)*1000, ParameterSet(qi,2), DummyEcho);  

    dual_Sat_val(qi) =calcMTsatThruLookupTablewithDummyV3( dual_GM_val(qi), ...
        1, T1, 1, M0, echoSpacing, ParameterSet(qi,6) + DummyEcho, ...
        ParameterSet(qi,3)*1000, ParameterSet(qi,2), DummyEcho);  

end

save( strcat(DATADIR,'Boost_single_Sat_val.mat'), 'single_Sat_val')
save( strcat(DATADIR,'Boost_dual_Sat_val.mat'),   'dual_Sat_val')






%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECTION NEEDS FILLED IN %%%%%%%%%%%%%%%%%%%%%%%
% Finish by calculating MTR, ihMTR and ihMTsat values. Combine into one matrix and save.


ihMT_sat=  dual_Sat_val - single_Sat_val ; 
sMTR  = (ref_GM_val - single_GM_val)./ref_GM_val;
dMTR  = (ref_GM_val - dual_GM_val)./ref_GM_val;
ihMTR = (single_GM_val - dual_GM_val)./ ref_GM_val;  


simResults = [single_GM_val, dual_GM_val, ref_GM_val, sMTR, dMTR, ihMTR,...
     single_Sat_val,dual_Sat_val, ihMT_sat];


save( strcat(DATADIR,'Boost_simResults_calcMetrics.mat'), 'simResults')

disp('Results table saved, script is complete!' )









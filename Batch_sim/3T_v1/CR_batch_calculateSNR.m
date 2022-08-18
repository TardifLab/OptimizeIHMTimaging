% Run this script after batch simulations, and before running the
% plotSimResultsIHMT_figures_SNR.m script. Running this script isn't
% necessary, but it should speed things up.

cl = parcluster('local');
cl.NumWorkers = 16;
saveProfile(cl);


%% I Will want to know this later...

% Contour lines drawn every SNR value of 1.

%% The goal of this is to generate plots 
% that show how: 1. contrast changes with different sequence parameters
% 2. how contrast efficiency changes with different sequence parameters. 

DATADIR = '/Directory/Where/Sim/Results/Are/Saved/';

%% Load not boosted values (conventional), and sequence parameters
convSim = load (strcat(DATADIR,'outputsConventional/nonBoost_simResults_calcMetrics.mat') );
convSim = convSim.simResults;

convParam = load (strcat(DATADIR,'outputsConventional/ParameterSet_3T_batch_conventional.mat') );
convParam = convParam.ParameterSet;

%% Load boosted , and sequence parameters
% and some 0 vectors and move the last two columns over
boostSim = load (strcat(DATADIR,'outputsBoost/Boost_simResults_calcMetrics.mat') );
boostSim = boostSim.simResults;

boostParam = load (strcat(DATADIR,'outputsBoost/ParameterSet_3T_batch_boost.mat') );
boostParam = boostParam.ParameterSet;


%% Load B1 values, they will be the same for single and dual, just load one
cB1 = load (strcat(DATADIR,'outputsConventional/nonBoost_dual_B1_val.mat') );
cB1 = cB1.dual_B1_val;

bB1 = load (strcat(DATADIR,'outputsBoost/Boost_dual_B1_val.mat') );
bB1 = bB1.dual_B1_val;

%% Concatenate

convParam  = [convParam cB1];
boostParam = [boostParam bB1];

simResults  = convSim;
simResultsB = boostSim;

clear cB1 bB1 convSim boostSim

% % for plotly
% I am getting a few with negative ihMT, it would be worth looking into
% those... Until then, don't do this. Looks like just poor saturation and
% excitation combos: 
% temp = convParam(simResults(:,9) <0,:);
% temp1 = simResults(simResults(:,9) <0,:);

% simResults  = abs(simResults);
% simResultsB = abs(simResultsB);

%% Prep a few extra values: 

% remove rows where SNR is too low and/or negative. 
s_thres = 0.035;
convParam(  simResults(:,1)  < s_thres,:) = [];
simResults( simResults(:,1)  < s_thres,:) = [];
boostParam( simResultsB(:,1) < s_thres,:) = [];
simResultsB(simResultsB(:,1) < s_thres,:) = [];

% Add a column that includes acquisition time.
% Setup image parameters for timing
numLines = 216*192; % matrix size % Mark Nelson Protocol
GRAPPA = 2;
refLines = 32;
PF = 1;
ellipticalScanning = 1; % set to 0 or 1
timeAcceptedMins = 4.5; % accepted time PER IMAGE

if ellipticalScanning
    numLines = round(((numLines/GRAPPA)+refLines) *PF * 0.78);
else
    numLines = ((numLines/GRAPPA)+refLines) *PF;
end

repit = numLines ./ (convParam(:,6));
timeAcq = repit .* convParam(:,3);
simResults(:,10) = timeAcq;

repit = numLines ./ (boostParam(:,6));
timeAcq = repit .* boostParam(:,3);
simResultsB(:,10) = timeAcq;

%% Compute SNR's if you haven't done previous (takes a while!)
parpool(str2num(getenv('NSLOTS')));

tic
if isfile(strcat(DATADIR,'outputsConventional/simResults_3T_batch_withCalc.mat'))
     
    simResults = load(strcat(DATADIR,'outputsConventional/simResults_3T_batch_withCalc.mat'));
    simResults = simResults.simResults;

else
     % File does not exist - do computation
    simResults = CR_add_ihMTR_SNR(simResults);
    disp('Finished Conventional MTR, moving to MTsat')
    qdewl
    simResults = CR_add_ihMTsat_SNR(simResults, convParam);
    
    save(strcat(DATADIR,'outputsConventional/simResults_3T_batch_withCalc.mat'),'simResults' );
end
toc

disp('Finished Conventional MT, moving to boosted')

% Boosted style
tic
if isfile(strcat(DATADIR,'outputsBoost/simResults_3T_batch_boost_withCalc.mat'))
     
    simResultsB = load(strcat(DATADIR,'outputsBoost/simResults_3T_batch_boost_withCalc.mat'));
    simResultsB = simResultsB.simResultsB;
else
     % File does not exist - do computation
    simResultsB = CR_add_ihMTR_SNR(simResultsB);
    disp('Finished Boosted MTR, moving to MTsat')
    
    simResultsB = CR_add_ihMTsat_SNR_boost(simResultsB, boostParam);
    
    save(strcat(DATADIR,'outputsBoost/simResults_3T_batch_boost_withCalc.mat'),'simResultsB' );
end
toc

disp('Finished Boosted MT, script complete.')




















%% The goal of this is to generate plots 
% that show how: 1. contrast changes with different sequence parameters
% 2. how contrast efficiency changes with different sequence parameters. 

DATADIR = '/Directory/Where/Sim/Results/Are/Saved/';
saveImgDir = '/Directory/To/Save/Output/Figures/';

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

clear cB1 bB1

simResults  = convSim;
simResultsB = boostSim;

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
tic
if isfile(strcat(DATADIR,'outputsConventional/simResults_3T_batch_withCalc.mat'))
     
    simResults = load(strcat(DATADIR,'outputsConventional/simResults_3T_batch_withCalc.mat'));
    simResults = simResults.simResults;

else
     % File does not exist - do computation
    simResults = CR_add_ihMTR_SNR(simResults);
    simResults = CR_add_ihMTsat_SNR(simResults, convParam);
    
    save(strcat(DATADIR,'outputsConventional/simResults_3T_batch_withCalc.mat'),'simResults' );
end
toc

% Boosted style
tic
if isfile(strcat(DATADIR,'outputsBoost/simResults_3T_batch_boost_withCalc.mat'))
     
    simResultsB = load(strcat(DATADIR,'outputsBoost/simResults_3T_batch_boost_withCalc.mat'));
    simResultsB = simResultsB.simResultsB;
else
     % File does not exist - do computation
    simResultsB = CR_add_ihMTR_SNR(simResultsB);
    simResultsB = CR_add_ihMTsat_SNR_boost(simResultsB, boostParam);
    
    save(strcat(DATADIR,'outputsBoost/simResults_3T_batch_boost_withCalc.mat'),'simResultsB' );
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% Export CSV for plotly
% convParamPlotly = array2table(convParam, 'VariableNames',{'delta', 'flipAngle', ...
%     'TR', 'numSatPulse', 'pulseDur', 'numExc', 'satTrainPerBoost', 'TR_MT','b1'});
% 
% boostParamPlotly = array2table(boostParam, 'VariableNames',{'delta', 'flipAngle', ...
%     'TR', 'numSatPulse', 'pulseDur', 'numExc', 'satTrainPerBoost', 'TR_MT','b1'});
% 
% simResultsPlotly = array2table(simResults, 'VariableNames',{ ...
%     'singleSig', 'dualSig', 'referenceSignal', ...
%     'singleMTR','dualMTR','ihMTR', ...
%     'MTsatSingle', 'MTsatDual', 'ihMTsat', ...
%     'acqTime', 'ihMTR_SNR','ihMTR_SNR_eff','ihMTsatSNR','ihMTsatSNReff'});
% 
% simResultsBPlotly = array2table(simResultsB, 'VariableNames',{ ...
%     'singleSig', 'dualSig', 'referenceSignal', ...
%     'singleMTR','dualMTR','ihMTR', ...
%     'MTsatSingle', 'MTsatDual', 'ihMTsat', ...
%     'acqTime', 'ihMTR_SNR','ihMTR_SNR_eff','ihMTsatSNR','ihMTsatSNReff'});
% 
% writetable(convParamPlotly, strcat(DATADIR,'outputsSNR/Plotly_convParams.csv')) 
% writetable(boostParamPlotly, strcat(DATADIR,'outputsSNR/Plotly_boostParams.csv')) 
% writetable(simResultsPlotly, strcat(DATADIR,'outputsSNR/Plotly_convSimResults.csv')) 
% writetable(simResultsBPlotly, strcat(DATADIR,'outputsSNR/Plotly_boostSimResults.csv')) 
% 
% clear convParamPlotly boostParamPlotly simResultsPlotly simResultsBPlotly



%% Column values in simResults
% single_GM_val, dual_GM_val, ref_GM_val, (1-3)
% sMTR, dMTR, ihMTR,...     (4-6)
% single_Sat_val,dual_Sat_val, ihMT_sat]; (7-9)
% Time, 'ihMTR_SNR','ihMTR_SNR_eff','ihMTsatSNR','ihMTsatSNReff' (10-14)

% convParam boostParam ->
% Offset_freq(a), flipAngle(b), TR(c), numSatPulse(d),...
% pulseDur(e), numExc(f), satTrainPerBoost(g), TR_MT(h), b1; 


%% Remove lines based on time
TF2 = simResults(:,10) > timeAcceptedMins*60; 
simResults(TF2,:) = [];
convParam(TF2,:) = [];

TF2 = simResultsB(:,10) > timeAcceptedMins*60; 
simResultsB(TF2,:) = [];
boostParam(TF2,:) = [];

%% Issues with some negative MTR values, remove those:
TF2 = (simResults(:,4) < 0) | (simResults(:,3) < 0) ; 
simResults(TF2,:) = [];
convParam(TF2,:) = [];

TF2 =  (simResultsB(:,4) < 0) | (simResultsB(:,3) < 0) ; 
simResultsB(TF2,:) = [];
boostParam(TF2,:) = [];

%% Remove lines based on signal
s_thres = 0.04;
convParam(  simResults(:,1)  < s_thres,:) = [];
simResults( simResults(:,1)  < s_thres,:) = [];
boostParam( simResultsB(:,1) < s_thres,:) = [];
simResultsB(simResultsB(:,1) < s_thres,:) = [];

%% Quick display of best sorted protocols:
%% Conventional
% sort by ihMTR SNR
% temp = simResults(top10Effidx,:);
[~, top10Effidx] = sort(simResults(:,11), 'descend'); % sort by ihMTsat SNR efficiency
Top10sorted = convParam ( top10Effidx,:);
Top10sorted = array2table(Top10sorted, 'VariableNames',{'delta', 'flipAngle', ...
    'TR', 'numSatPulse', 'pulseDur', 'numExc', 'satTrainPerBoost', 'TR_MT','b1'});

% REDO above but for ihMTsat SNR
[~, top10Effidx] = sort(simResults(:,13), 'descend'); % sort by ihMTR SNR efficiency
Top10sortedSat = convParam ( top10Effidx(1:20),:);
Top10sortedSat = array2table(Top10sortedSat, 'VariableNames',{'delta', 'flipAngle', ...
    'TR', 'numSatPulse', 'pulseDur', 'numExc', 'satTrainPerBoost', 'TR_MT','b1'});


%% Boosted
% sort by ihMTR SNR
temp = simResultsB(top10Effidx,:);
[~, top10Effidx] = sort(simResultsB(:,11), 'descend'); % sort by ihMTsat SNR efficiency
Top10sorted_b = boostParam ( top10Effidx,:);
Top10sorted_b = array2table(Top10sorted_b, 'VariableNames',{'delta', 'flipAngle', ...
    'TR', 'numSatPulse', 'pulseDur', 'numExc', 'satTrainPerBoost', 'TR_MT','b1'});

% REDO above but for ihMTsat SNR
[~, top10Effidx] = sort(simResultsB(:,13), 'descend'); % sort by ihMTR SNR efficiency
Top10sortedSat_b = boostParam ( top10Effidx(1:20),:);
Top10sortedSat_b = array2table(Top10sortedSat_b, 'VariableNames',{'delta', 'flipAngle', ...
    'TR', 'numSatPulse', 'pulseDur', 'numExc', 'satTrainPerBoost', 'TR_MT','b1'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% We are interested in speeding it up. 
% Generate a 3D plot, TRvs numExcitations vs ihMTsat contrast efficiency.

[TRm, numExcm, SNReff, SNRsateff,SNRabs, SNRsatabs] = ...
    ihMTsim_SNRvectorGen4Figures( convParam, [3,6], 11:14, simResults);


[TRmB, numExcmB, SNReffB, SNRsateffB,SNRabsB, SNRsatabsB] = ...
    ihMTsim_SNRvectorGen4Figures( boostParam, [3,6], 11:14, simResultsB);


%%
absLowS = 0; absHighS = max(SNRabs,[],'all') + 0.1*max(SNRabs,[],'all');
absLowBS = 0; absHighBS = max(SNRabsB,[],'all') + 0.1*max(SNRabsB,[],'all');

absLowSatS = 0; absHighSatS = max(SNRsatabs,[],'all') + 0.1*max(SNRsatabs,[],'all');
absLowSatBS = 0; absHighSatBS = max(SNRsatabsB,[],'all') + 0.1*max(SNRsatabsB,[],'all');

maxPlot = max([absHighS, absHighBS, absHighSatS,absHighSatBS]);

%% -> Separate for this one.
ihMTsimResultPlotting_v2 (TRm, numExcm, SNRabs,...
    'Turbofactor',    'TR (sec)',    'ihMTR SNR',...
    absLowS, maxPlot, 'surf' , 0)
saveas(gcf,strcat(saveImgDir,'ExcitationsVsTRvsihMTRsnr.png'))

ihMTsimResultPlotting_v2 (TRm, numExcm, SNRsatabs,...
    'Turbofactor',    'TR (sec)',    'ihMT_{sat} SNR',...
    absLowSatS, maxPlot, 'surf' , 0)
saveas(gcf,strcat(saveImgDir,'ExcitationsVsTRvsihMTsatsnr.png'))
    
%% Boosted 

ihMTsimResultPlotting_v2 (TRmB, numExcmB, SNRabsB,...
    'Turbofactor',    'TR (sec)',    'ihMTR',...
    absLowBS, maxPlot, 'surf' , 0)
saveas(gcf,strcat(saveImgDir,'ExcitationsVsTRvsihMTRBsnr.png'))


ihMTsimResultPlotting_v2 (TRmB, numExcmB, SNRsatabsB,...
    'Turbofactor',    'TR (sec)',    'ihMT_{sat} SNR',...
    absLowSatBS, maxPlot, 'surf' , 0)
saveas(gcf,strcat(saveImgDir,'ExcitationsVsTRvsihMTsatBsnr.png'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare TR to flipAngle


[TRm, flipAnglem, SNReff, SNRsateff, SNRabs, SNRsatabs] = ...
    ihMTsim_SNRvectorGen4Figures( convParam, [3,2], 11:14, simResults);


[TRmB, flipAnglemB, SNReffB, SNRsateffB, SNRabsB, SNRsatabsB] = ...
    ihMTsim_SNRvectorGen4Figures( boostParam, [3,2], 11:14, simResultsB);


ihMTsimResultPlotting_v2 (TRm, flipAnglem, SNRabs,...
    'Flip Angle (deg)',    'TR (sec)',    'ihMTR SNR',...
    absLowS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'flipAngleVsTRvsihMTRsnr.png'))


ihMTsimResultPlotting_v2 (TRm, flipAnglem, SNRsatabs,...
    'Flip Angle (deg)',    'TR (sec)',    'ihMT_{sat} SNR',...
    absLowSatS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'flipAngleVsTRvsihMTsatsnr.png'))


%% Boosted

ihMTsimResultPlotting_v2 (TRmB, flipAnglemB, SNRabsB,...
    'Flip Angle (deg)',    'TR (sec)',    'ihMTR SNR',...
    absLowBS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'flipAngleVsTRvsihMTRBsnr.png'))


ihMTsimResultPlotting_v2 (TRmB, flipAnglemB, SNRsatabsB,...
    'Flip Angle (deg)',    'TR (sec)',    'ihMT_{sat} SNR',...
    absLowSatBS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'flipAngleVsTRvsihMTsatBsnr.png'))




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare TR to number of SatPulses

[TRm, numSatPulm, SNReff, SNRsateff, SNRabs, SNRsatabs] = ...
    ihMTsim_SNRvectorGen4Figures( convParam, [3,4], 11:14, simResults);


[TRmB, numSatPulmB, SNReffB, SNRsateffB, SNRabsB, SNRsatabsB] = ...
    ihMTsim_SNRvectorGen4Figures( boostParam, [3,4], 11:14, simResultsB);


ihMTsimResultPlotting_v2 (TRm, numSatPulm, SNRabs,...
    'N_{sat}',...
    'TR (sec)',...
    'ihMTR SNR',...
    absLowS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'numSatPulVsTRvsihMTRsnr.png'))


ihMTsimResultPlotting_v2 (TRm, numSatPulm, SNRsatabs,...
    'N_{sat}',    'TR (sec)',    'ihMT_{sat} SNR',...
    absLowSatS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'numSatPulVsTRvsihMTsatsnr.png'))



%% BOOSTED

ihMTsimResultPlotting_v2 (TRmB, numSatPulmB, SNRabsB,...
    'N_{sat}',...
    'TR (sec)',...
    'ihMTR SNR',...
    absLowBS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'numSatPulVsTRvsihMTRsnrB.png'))


ihMTsimResultPlotting_v2 (TRmB, numSatPulmB, SNRsatabsB,...
    'N_{sat}',    'TR (sec)',    'ihMT_{sat} SNR',...
    absLowSatBS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'numSatPulVsTRvsihMTsatBsnr.png'))




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Compare numExcitation to number of SatPulses 

[numSatPulm, numExcm, SNReff, SNRsateff, SNRabs, SNRsatabs] = ...
    ihMTsim_SNRvectorGen4Figures( convParam, [4,6], 11:14, simResults);


[numSatPulmB, numExcmB, SNReffB, SNRsateffB, SNRabsB, SNRsatabsB] = ...
    ihMTsim_SNRvectorGen4Figures( boostParam, [4,6], 11:14, simResultsB);



ihMTsimResultPlotting_v2 (numSatPulm, numExcm, SNRabs,...
    'Turbofactor',...
    'N_{sat}',...
    'ihMTR SNR',...
    absLowS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'ExcitationsVsnumSatPulvsihMTRsnr.png'))


ihMTsimResultPlotting_v2 (numSatPulm, numExcm, SNRsatabs,...
    'Turbofactor',...
    'N_{sat}',...
    'ihMT_{sat} SNR',...
    absLowSatS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'ExcitationsVsnumSatPulvsihMTsatsnr.png'))



%% BOOST  
   
ihMTsimResultPlotting_v2 (numSatPulmB, numExcmB, SNRabsB,...
    'Turbofactor',...
    'N_{sat}',...
    'ihMTR SNR',...
    absLowBS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'ExcitationsVsnumSatPulvsihMTRBsnr.png'))


ihMTsimResultPlotting_v2 (numSatPulmB, numExcmB, SNRsatabsB,...
    'Turbofactor',...
    'N_{sat}',...
    'ihMT_{sat} SNR',...
    absLowSatBS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'ExcitationsVsnumSatPulvsihMTsatBsnr.png'))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
 %% Compare PulseDur to numSatPulse    

 [numSatPulm, pulseDurm, SNReff, SNRsateff, SNRabs, SNRsatabs] = ...
    ihMTsim_SNRvectorGen4Figures( convParam, [4,5], 11:14, simResults);


[numSatPulmB, pulseDurmB, SNReffB, SNRsateffB, SNRabsB, SNRsatabsB] = ...
    ihMTsim_SNRvectorGen4Figures( boostParam, [4,5], 11:14, simResultsB);
 
 
 
ihMTsimResultPlotting_v2 (numSatPulm, pulseDurm*1000, SNRabs,...
    't_{pulse} (ms)',...
   'N_{sat}',...
    'ihMTR SNR',...
    absLowS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'SatPulseDurVsnumSatPulvsihMTRsnr.png'))


ihMTsimResultPlotting_v2 (numSatPulm, pulseDurm*1000, SNRsatabs,...
   't_{pulse} (ms)',...
   'N_{sat}',...
   'ihMT_{sat} SNR',...
    absLowSatS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'SatPulseDurVsnumSatPulvsihMTsatsnr.png'))

%% BOOST
 
ihMTsimResultPlotting_v2 (numSatPulmB, pulseDurmB*1000, SNRabsB,...
    't_{pulse} (ms)',...
   'N_{sat}',...
    'ihMTR SNR',...
    absLowBS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'SatPulseDurVsnumSatPulvsihMTRBsnr.png'))


ihMTsimResultPlotting_v2 (numSatPulmB, pulseDurmB*1000, SNRsatabsB,...
   't_{pulse} (ms)',...
   'N_{sat}',...
   'ihMT_{sat} SNR',...
    absLowSatBS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'SatPulseDurVsnumSatPulvsihMTsatBsnr.png'))




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Compare satTrainPerBoost to TR_MT
% Boosted Protocol only!

 [TR_MTmB, satTrainPerBoostmB, SNReffB, SNRsateffB, SNRabsB, SNRsatabsB] = ...
    ihMTsim_SNRvectorGen4Figures( boostParam, [8,7], 11:14, simResultsB);


 
ihMTsimResultPlotting_v2 (TR_MTmB*1000, satTrainPerBoostmB, SNRabsB,...
    'N_{boost}', 'TR_{MT} (ms)',   'ihMTR SNR',...
    absLowBS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'TR_MTvsSatTrainPerBoostvsihMTRBsnr.png'))


ihMTsimResultPlotting_v2 (TR_MTmB*1000, satTrainPerBoostmB, SNRsatabsB,...
    'N_{boost}', 'TR_{MT} (ms)',    'ihMT_{sat} SNR',...
    absLowSatBS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'TR_MTvsSatTrainPerBoostvsihMTsatBsnr.png'))





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Compare N_sat to N_boost
% Boosted Protocol only!

 
 [numSatPulm, satTrainPerBoostmB, SNReffB, SNRsateffB, SNRabsB, SNRsatabsB] = ...
    ihMTsim_SNRvectorGen4Figures( boostParam, [4,7], 11:14, simResultsB);

 
ihMTsimResultPlotting_v2 (numSatPulm, satTrainPerBoostmB, SNRabsB,...
   'N_{boost}', 'N_{sat}',  'ihMTR SNR',...
    absLowBS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'NsatVsNboostVsihMTRBsnr.png'))


ihMTsimResultPlotting_v2 (numSatPulm, satTrainPerBoostmB, SNRsatabsB,...
    'N_{boost}', 'N_{sat}',    'ihMT_{sat} SNR',...
    absLowSatBS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'NsatVsNboostvsihMTsatBsnr.png'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% Compare N_sat to TR_MT
% Boosted Protocol only!

[TR_MTmB, numSatPulm, SNReffB, SNRsateffB, SNRabsB, SNRsatabsB] = ...
    ihMTsim_SNRvectorGen4Figures( boostParam, [8,4], 11:14, simResultsB);
 
ihMTsimResultPlotting_v2 (TR_MTmB*1000, numSatPulm, SNRabsB,...
   'N_{sat}','TR_{MT} (ms)',  'ihMTR SNR',...
    absLowBS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'TR_MTVsnumSatPulVsihMTRBsnr.png'))



ihMTsimResultPlotting_v2 (TR_MTmB*1000, numSatPulm, SNRsatabsB,...
    'N_{sat}','TR_{MT} (ms)',    'ihMT_{sat} SNR',...
    absLowSatBS, maxPlot, 'surf' , 1)
saveas(gcf,strcat(saveImgDir,'TR_MTVsnumSatPulvsihMTsatBsnr.png'))






                




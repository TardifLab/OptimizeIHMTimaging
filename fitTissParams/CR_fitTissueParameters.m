
%% Estimate tissue parameters for B1 correction simulations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baseDir = '/path/to/OptimizeIHMTimaging/';
saveImgDir = '/Directory/To/Save/Output/Figures/';


%% Load data
ref_kspace_dir = strcat(baseDir,'kspaceWeighting/Atlas_reference/');
load( strcat( ref_kspace_dir,'GM_seg_MNI_152_image.mat')) % keep these here so you know what files you are setting directory to
load( strcat( ref_kspace_dir,'GM_seg_MNI_152_kspace.mat'))

load( [baseDir,'fitTissParams/MTsat_vals2fit.mat'] )   
load( [baseDir,'fitTissParams/MTsat_vals2fit_Mar16.mat'] )   

%% Set up
fit_version = '1p0_FitParams'; % used as naming convention

SavDir = [SavDir,fit_version,'/'];
mkdir(SavDir)

%% Start including a change log:
logString= strcat('Notes for this fit');

fid = fopen(strcat(SavDir,'log.txt'),'wt');
fprintf(fid, logString);
fclose(fid);


R = linspace(15,50,4);
T2a = linspace(20e-3,90e-3,3);
T1D =  [5e-4 1e-3 5e-3];% Varma 2017 was 6ms
T2b = linspace(8e-6, 12e-6, 3);
%M0b =  linspace(0.0475, 0.06, 5);  % include this to give all parameters a chance...
M0b = [0.065, 0.07, 0.075]; 
R1b = [0.25, 0.75, 1];  % can do brief sims with this one at end... Has very little impact on its own.
B1rms = 0:2:18;

%simLength = length(R)* length(T2a)* length(T1D)* length(T2b)* length(M0b)* length(R1b);
simLength2 = length(R)* length(T2a)* length(T1D)* length(T2b)* length(M0b)* length(R1b);

%%
% First set a couple of initial params, then fill defaults, then set the
% rest. 
Params.B0 = 3;
Params.MTC = 1; % Magnetization Transfer Contrast
Params.TissueType = 'GM';
Params = DefaultCortexTissueParams(Params);
Params = CalcImagingParams(Params);


Params.b1 = 0; % microTesla
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
Params.PerfectSpoiling = 1;

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
[outputSamplingTable, ~, Params.Segments] = Step1_calculateKspaceSampling_v3 (Params);


%% Other 2 sequences:

Params2 = Params;
Params3 = Params;

Params2.b1 = 0; % microTesla
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
[outputSamplingTable2, ~, Params2.Segments] = Step1_calculateKspaceSampling_v3 (Params2);


Params3.b1 = 0; % microTesla
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
[outputSamplingTable3, ~, Params3.Segments] = Step1_calculateKspaceSampling_v3 (Params3);



%% Tissue parameters

%% Setup parameter structure
parametersSet2 = zeros( simLength2,  6);

idx = 1;
for a = 1:length(R)    
    for b = 1:length(T2a) 
        for c = 1:length(T1D)
            for d = 1:length(T2b)
                for e = 1:length(M0b)
                    for f = 1:length(R1b)
                        
                        parametersSet2(idx,:) = [R(a),...
                            T2a(b), T1D(c), T2b(d),...
                            M0b(e), R1b(f) ]; 

                        idx = idx+1;
                    end
                end
            end
        end
    end
end

save(strcat(SavDir,'parametersSet.mat'),'parametersSet2') % save this just incase I change it at some point...




%% Run simulations
TF1 = Params.TurboFactor;
TF2 = Params2.TurboFactor;
TF3 = Params3.TurboFactor;

Single_sig1 = zeros( simLength2, TF1, length(B1rms));
Single_sig2 = zeros( simLength2, TF2, length(B1rms));
Single_sig3 = zeros( simLength2, TF3, length(B1rms));
Dual_sig1 = zeros( simLength2, TF1,  length(B1rms));
Dual_sig2 = zeros( simLength2, TF2, length(B1rms));
Dual_sig3 = zeros( simLength2, TF3, length(B1rms));

% run the
parpool

tic % Took 90 hours to run

parfor qi = 1:simLength2

    t1 = parametersSet2(qi,1);
    t2 = parametersSet2(qi,2);
    t3 = parametersSet2(qi,3);
    t4 = parametersSet2(qi,4);
    t5 = parametersSet2(qi,5);
    t6 = parametersSet2(qi,6);

    Ssig1 = zeros(TF1, length(B1rms) );
    Ssig2 = zeros(TF2, length(B1rms) );
    Ssig3 = zeros(TF3, length(B1rms) );
    Dsig1 = zeros(TF1, length(B1rms) );
    Dsig2 = zeros(TF2, length(B1rms) );
    Dsig3 = zeros(TF3, length(B1rms) );

    for i = 1:10
         [Ssig1(:,i),~, ~]  = BlochSimFlashSequence_v2(Params,'freqPattern', 'single', 'b1', B1rms(i),...
             'R', t1, 'T2a', t2, 'T1D',t3, 'T2b', t4, 'M0b', t5, 'R1b', t6 ); 
                        
         [Ssig2(:,i),~, ~]  = BlochSimFlashSequence_v2(Params2,'freqPattern', 'single', 'b1', B1rms(i),...
             'R', t1, 'T2a', t2, 'T1D',t3, 'T2b', t4, 'M0b', t5, 'R1b', t6 );
         
         [Ssig3(:,i),~, ~]  = BlochSimFlashSequence_v2(Params3,'freqPattern', 'single', 'b1', B1rms(i),...
             'R', t1, 'T2a', t2, 'T1D',t3, 'T2b', t4, 'M0b', t5, 'R1b', t6 );
         
         [Dsig1(:,i),~, ~]  = BlochSimFlashSequence_v2(Params,'freqPattern', 'dualAlternate', 'b1', B1rms(i),...
             'R', t1, 'T2a', t2, 'T1D',t3, 'T2b', t4, 'M0b', t5, 'R1b', t6 );
              
         [Dsig2(:,i),~, ~]  = BlochSimFlashSequence_v2(Params2,'freqPattern', 'dualAlternate', 'b1', B1rms(i),...
             'R', t1, 'T2a', t2, 'T1D',t3, 'T2b', t4, 'M0b', t5, 'R1b', t6 );

         [Dsig3(:,i),~, ~]  = BlochSimFlashSequence_v2(Params3,'freqPattern', 'dualAlternate', 'b1', B1rms(i),...
             'R', t1, 'T2a', t2, 'T1D',t3, 'T2b', t4, 'M0b', t5, 'R1b', t6 ); 
    end
                        
    Single_sig1(qi,:,:) = Ssig1; 
    Single_sig2(qi,:,:) = Ssig2;
    Single_sig3(qi,:,:) = Ssig3; 
    Dual_sig1(qi,:,:)   = Dsig1;
    Dual_sig2(qi,:,:)   = Dsig2;
    Dual_sig3(qi,:,:)   = Dsig3;
                                               
    if rem(qi, 50) == 0

        qi/simLength2 *100 % print percent done
    end   
end
toc

save(strcat(SavDir,'Single_sig1.mat'),'Single_sig1')
save(strcat(SavDir,'Single_sig2.mat'),'Single_sig2')
save(strcat(SavDir,'Single_sig3.mat'),'Single_sig3')
save(strcat(SavDir,'Dual_sig1.mat'),'Dual_sig1')
save(strcat(SavDir,'Dual_sig2.mat'),'Dual_sig2')
save(strcat(SavDir,'Dual_sig3.mat'),'Dual_sig3')



%% Then from the excitation train values, determine the realized GM value.

Single_sig1_gm = zeros( simLength2, length(B1rms));
Single_sig2_gm = zeros( simLength2, length(B1rms));
Single_sig3_gm = zeros( simLength2, length(B1rms));
Dual_sig1_gm = zeros( simLength2,  length(B1rms));
Dual_sig2_gm = zeros( simLength2, length(B1rms));
Dual_sig3_gm = zeros( simLength2,  length(B1rms));

for i = 1:simLength2
    %parfor j = 1:length(B1rms)
    for j = 1:length(B1rms)

        Single_sig1_gm(i,j) = CR_generate_BSF_scaling_v1( squeeze(Single_sig1(i,:,j)), Params, outputSamplingTable, gm_m, fft_gm_m) ;  
        Single_sig2_gm(i,j) = CR_generate_BSF_scaling_v1(squeeze(Single_sig2(i,:,j)), Params2, outputSamplingTable2, gm_m, fft_gm_m) ;
        Single_sig3_gm(i,j) = CR_generate_BSF_scaling_v1(squeeze(Single_sig3(i,:,j)), Params3, outputSamplingTable3, gm_m, fft_gm_m) ;
        Dual_sig1_gm(i,j) = CR_generate_BSF_scaling_v1(squeeze(Dual_sig1(i,:,j)), Params, outputSamplingTable, gm_m, fft_gm_m) ;
        Dual_sig2_gm(i,j) = CR_generate_BSF_scaling_v1(squeeze(Dual_sig2(i,:,j)), Params2, outputSamplingTable2, gm_m, fft_gm_m);
        Dual_sig3_gm(i,j) = CR_generate_BSF_scaling_v1(squeeze(Dual_sig3(i,:,j)), Params3, outputSamplingTable3, gm_m, fft_gm_m) ;  

    end
end


save(strcat(SavDir,'Single_sig1_gm.mat'),'Single_sig1_gm')
save(strcat(SavDir,'Single_sig2_gm.mat'),'Single_sig2_gm')
save(strcat(SavDir,'Single_sig3_gm.mat'),'Single_sig3_gm')
save(strcat(SavDir,'Dual_sig1_gm.mat'),'Dual_sig1_gm')
save(strcat(SavDir,'Dual_sig2_gm.mat'),'Dual_sig2_gm')
save(strcat(SavDir,'Dual_sig3_gm.mat'),'Dual_sig3_gm')


%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculate MTsat on the whole matrix of signal values.
T1obs = ones(size(Dual_sig3_gm)) .* 1.4.*1000;
M0_app_v = ones(size(Dual_sig3_gm)) ;


MTsat_sim_Single1 = calcMTsatThruLookupTablewithDummyV3( Single_sig1_gm, [], T1obs, [], M0_app_v, Params.echoSpacing * 1000,  Params.numExcitation,  Params.TR * 1000,  Params.flipAngle,   Params.DummyEcho);
MTsat_sim_Single2 = calcMTsatThruLookupTablewithDummyV3( Single_sig2_gm, [], T1obs, [], M0_app_v, Params2.echoSpacing * 1000, Params2.numExcitation, Params2.TR * 1000, Params2.flipAngle, Params2.DummyEcho);
MTsat_sim_Single3 = calcMTsatThruLookupTablewithDummyV3( Single_sig3_gm, [], T1obs, [], M0_app_v, Params3.echoSpacing * 1000, Params3.numExcitation, Params3.TR * 1000, Params3.flipAngle, Params3.DummyEcho);
MTsat_sim_Dual1   = calcMTsatThruLookupTablewithDummyV3( Dual_sig1_gm,   [], T1obs, [], M0_app_v, Params.echoSpacing * 1000,  Params.numExcitation,  Params.TR * 1000,  Params.flipAngle,   Params.DummyEcho);
MTsat_sim_Dual2   = calcMTsatThruLookupTablewithDummyV3( Dual_sig2_gm,   [], T1obs, [], M0_app_v, Params2.echoSpacing * 1000, Params2.numExcitation, Params2.TR * 1000, Params2.flipAngle, Params2.DummyEcho);
MTsat_sim_Dual3   = calcMTsatThruLookupTablewithDummyV3( Dual_sig3_gm,   [], T1obs, [], M0_app_v, Params3.echoSpacing * 1000, Params3.numExcitation, Params3.TR * 1000, Params3.flipAngle, Params3.DummyEcho);




%% Remove row zeros to save time:
maxRowS = max(MTsat_sim_Dual1,[],2);

MTsat_sim_Single1(maxRowS == 0,:) = [];
MTsat_sim_Single2(maxRowS == 0,:) = [];
MTsat_sim_Single3(maxRowS == 0,:) = [];
MTsat_sim_Dual1(maxRowS == 0,:) = [];
MTsat_sim_Dual2(maxRowS == 0,:) = [];
MTsat_sim_Dual3(maxRowS == 0,:) = [];
parametersSet2(maxRowS == 0,:) = [];

%% Fit the data
fit_degree = 7;
Single_c1 = zeros( length(MTsat_sim_Single1),  fit_degree+1);
Single_c2 = zeros( length(MTsat_sim_Single1),  fit_degree+1);
Single_c3 = zeros( length(MTsat_sim_Single1),  fit_degree+1);
Dual_c1 = zeros( length(MTsat_sim_Dual1), fit_degree+1);
Dual_c2 = zeros( length(MTsat_sim_Dual1), fit_degree+1);
Dual_c3 = zeros( length(MTsat_sim_Dual1), fit_degree+1);

tic % super fast, few seconds  
for i = 1:length(MTsat_sim_Dual1)
                    
     % With B1's simulated ->Fit polynomial: 
     Single_c1(i,:) = polyfit(B1rms, MTsat_sim_Single1(i,:), fit_degree);
     Single_c2(i,:) = polyfit(B1rms, MTsat_sim_Single2(i,:), fit_degree);
     Single_c3(i,:) = polyfit(B1rms, MTsat_sim_Single3(i,:), fit_degree);
     Dual_c1(i,:)   = polyfit(B1rms, MTsat_sim_Dual1(i,:), fit_degree);
     Dual_c2(i,:)   = polyfit(B1rms, MTsat_sim_Dual2(i,:), fit_degree);
     Dual_c3(i,:)   = polyfit(B1rms, MTsat_sim_Dual3(i,:), fit_degree);
     
     
end
toc

% % Can do a quick check if you want!
% i = 1;
% x1 = linspace(0,18,100);
% y1 = polyval(Single_c1(i,:),x1);
% figure
% plot(B1rms, MTsat_sim_Single1(i,:),'o')
% hold on
% plot(x1,y1)
% hold off

                         
 %% With fitted data, calculate residuals: https://www.mathworks.com/help/matlab/ref/polyfit.html
 % To see how good the fit is, evaluate the polynomial at the data points and generate a table showing the data, fit, and error.      


% stack matrices:
% Sort B1 , then sat values, then smooth them
[~, sortidx_b1] = sort(  exportMat(10,:) , 'ascend');

% Sort the whole matrix:
mat_sort = exportMat( :, sortidx_b1);

% Smooth matrix
mat_ss = smoothdata( mat_sort,2, 'movmedian', 10);

% Make sure you have B1rms values! Multiply B1 map by the B1rms of the sequence 
b1_1 =11.4* mat_ss(10,:) ; 
b1_2 = 13.3 * mat_ss(10,:) ; 
b1_3 = 11.6*mat_ss(10,:) ; 

%% Quick plot of each to see that the data looks OK
% figure; heatscatter(b1_1', mat_ss(1,:)'); %ylim([0 0.04]); xlim([0 15])
% figure; heatscatter(b1_3', mat_ss(6,:)'); %ylim([0 0.04]); xlim([0 15])


%% Calculate Standardized Residuals
% this can take a few minutes of run time.
std_resid_d1 = CR_calc_std_residuals( b1_1, mat_ss(1,:) , Dual_c1);
std_resid_d2 = CR_calc_std_residuals( b1_2, mat_ss(2,:) , Dual_c2);
std_resid_d3 = CR_calc_std_residuals( b1_3, mat_ss(3,:) , Dual_c3);

std_resid_s1 = CR_calc_std_residuals( b1_1, mat_ss(4,:) , Single_c1);
std_resid_s2 = CR_calc_std_residuals( b1_2, mat_ss(5,:) , Single_c2);
std_resid_s3 = CR_calc_std_residuals( b1_3, mat_ss(6,:) , Single_c3);

standardized_residuals = [std_resid_d1, std_resid_d2, std_resid_d3, std_resid_s1, std_resid_s2, std_resid_s3];

 save(strcat(SavDir,'standardized_residuals.mat'),'standardized_residuals')        
 
 
 %%  add the columns
 errorScore = sum( abs(standardized_residuals) , 2 ); % can have negatives, so combine abs of each
 
 % Sort based on this, then store in table :) 
 
% for best protocol, take top 10 most efficient protocols. Then sort by
% absolute ihMT

[~, sortidx] = sort(errorScore, 'ascend');
Top50sorted = parametersSet2 ( sortidx(1:end),:);
Top50Errors = errorScore ( sortidx(1:end),:);
Top50sortedTable = array2table(Top50sorted, 'VariableNames',{'R', 'T2a', 'T1D', 'T2b', 'M0b','R1b'});

save(strcat(SavDir,'Top50sortedTable.mat'),'Top50sortedTable')      


str = ['R = ',num2str(Top50sorted(1,1)),', T2a = ',num2str(Top50sorted(1,2)),...
    ', T1D = ',num2str(Top50sorted(1,3)), ', T2b =',num2str(Top50sorted(1,4)),...
    ', M0b = ',num2str(Top50sorted(1,5)),', R1b = ',num2str(Top50sorted(1,6))];
                         
% Check what the top one looks like with the data!      
  % Check first few, number 1 looked off here, 2 was better
  
SortIndex = sortidx(1); % select the sorted line you want
x1_line = linspace(0,18,100);
y1 = polyval( Dual_c1(SortIndex,:), x1_line);
y2 = polyval( Dual_c2(SortIndex,:), x1_line);
y3 = polyval( Dual_c3(SortIndex,:), x1_line);



figure;
subplot(1,2,1)
heatscatter(b1_1', mat_ss(1,:)' ); 
hold on
heatscatter(b1_2', mat_ss(2,:)' ); 
heatscatter(b1_3', mat_ss(3,:)' ); 
hold on
plot(x1_line,y1,'LineWidth',3); plot(x1_line,y2,'LineWidth',3); plot(x1_line,y3,'LineWidth',3)
xlim([0 18]) ; % ylim([0 0.04]) ;
title('Dual');
ylabel('MT_{sat}); xlabel('B_{1rms}');
colorbar off
ax = gca; ax.FontSize = 20; 
hold off


y1 = polyval( Single_c1( SortIndex,:), x1_line);
y2 = polyval( Single_c2( SortIndex,:), x1_line);
y3 = polyval( Single_c3( SortIndex,:), x1_line);

subplot(1,2,2)
heatscatter(b1_1', mat_ss(4,:)' ); 
hold on
heatscatter(b1_2', mat_ss(5,:)' ); 
heatscatter(b1_3', mat_ss(6,:)' ); 
hold on
plot(x1_line,y1,'LineWidth',3); plot(x1_line,y2,'LineWidth',3); plot(x1_line,y3,'LineWidth',3)
xlim([0 18]) ; % ylim([0 0.04]) ;
title('Single')
ylabel('MT_{sat}); xlabel('B_{1rms}');
ax = gca; ax.FontSize = 20; 
colorbar off
hold off                     
  set(gcf,'position',[10,400,1200,400])                       
     sgtitle(str)             
     
     

   saveas(gcf,strcat(SavDir,'initial_best_parameters_fit2Optimal.png'))  
     
     
     
   
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %% Redo fit with second set of data.

% stack matrices:
% Sort B1 , then sat values, then smooth them
[~, sortidx_b1] = sort(  exportMat2(10,:) , 'ascend');

% Sort the whole matrix:
mat_sort = exportMat2( :, sortidx_b1);

% Smooth matrix
mat_ss = smoothdata( mat_sort,2, 'movmedian', 10);

% Make sure you have B1rms values! Multiply B1 map by the B1rms of the sequence 
b1_1 =11.4* mat_ss(10,:) ; 
b1_2 = 13.3 * mat_ss(10,:) ; 
b1_3 = 11.6*mat_ss(10,:) ; 

%% Quick plot of each to see that the data looks OK
% figure; heatscatter(b1_1', mat_ss(1,:)'); %ylim([0 0.04]); xlim([0 15])
% figure; heatscatter(b1_3', mat_ss(6,:)'); %ylim([0 0.04]); xlim([0 15])


%% Calculate Standardized Residuals
% this can take a few minutes of run time.
std_resid_d1 = CR_calc_std_residuals( b1_1, mat_ss(1,:) , Dual_c1);
std_resid_d2 = CR_calc_std_residuals( b1_2, mat_ss(2,:) , Dual_c2);
std_resid_d3 = CR_calc_std_residuals( b1_3, mat_ss(3,:) , Dual_c3);

std_resid_s1 = CR_calc_std_residuals( b1_1, mat_ss(4,:) , Single_c1);
std_resid_s2 = CR_calc_std_residuals( b1_2, mat_ss(5,:) , Single_c2);
std_resid_s3 = CR_calc_std_residuals( b1_3, mat_ss(6,:) , Single_c3);

standardized_residuals2 = [std_resid_d1, std_resid_d2, std_resid_d3, std_resid_s1, std_resid_s2, std_resid_s3];

save(strcat(SavDir,'standardized_residuals2.mat'),'standardized_residuals2')        
 


 %%  add the columns
 errorScore = sum( abs(standardized_residuals2) , 2 ); % can have negatives, so combine abs of each
 
 % Sort based on this, then store in table :) 
 
% for best protocol, take top 10 most efficient protocols. Then sort by
% absolute ihMT

[temp, sortidx2] = sort(errorScore, 'ascend');
Top50sorted2 = parametersSet2 ( sortidx2(1:end),:);
Top50Errors2 = errorScore ( sortidx2(1:end),:);
Top50sortedTable2 = array2table(Top50sorted2, 'VariableNames',{'R', 'T2a', 'T1D', 'T2b', 'M0b','R1b'});

save(strcat(SavDir,'Top50sortedTable2.mat'),'Top50sortedTable2')   


str = ['R = ',num2str(Top50sorted2(1,1)),', T2a = ',num2str(Top50sorted2(1,2)),...
    ', T1D = ',num2str(Top50sorted2(1,3)), ', T2b =',num2str(Top50sorted2(1,4)),...
    ', M0b = ',num2str(Top50sorted2(1,5)),', R1b = ',num2str(Top50sorted2(1,6))];
                         
% Check what the top one looks like with the data!      
  % Check first few, number 1 looked off here, 2 was better
  
SortIndex = sortidx2(1); % select the sorted line you want
x1_line = linspace(0,18,100);
y1 = polyval( Dual_c1(SortIndex,:), x1_line);
y2 = polyval( Dual_c2(SortIndex,:), x1_line);
y3 = polyval( Dual_c3(SortIndex,:), x1_line);



figure;
subplot(1,2,1)
heatscatter(b1_1', mat_ss(1,:)' ); 
hold on
heatscatter(b1_2', mat_ss(2,:)' ); 
heatscatter(b1_3', mat_ss(3,:)' ); 
hold on
plot(x1_line,y1,'LineWidth',3); plot(x1_line,y2,'LineWidth',3); plot(x1_line,y3,'LineWidth',3)
xlim([0 18]) ; % ylim([0 0.04]) ;
title('Dual');
ylabel('MT_{sat}); xlabel('B_{1rms}');
colorbar off
ax = gca; ax.FontSize = 20; 
hold off


y1 = polyval( Single_c1( SortIndex,:), x1_line);
y2 = polyval( Single_c2( SortIndex,:), x1_line);
y3 = polyval( Single_c3( SortIndex,:), x1_line);

subplot(1,2,2)
heatscatter(b1_1', mat_ss(4,:)' ); 
hold on
heatscatter(b1_2', mat_ss(5,:)' ); 
heatscatter(b1_3', mat_ss(6,:)' ); 
hold on
plot(x1_line,y1,'LineWidth',3); plot(x1_line,y2,'LineWidth',3); plot(x1_line,y3,'LineWidth',3)
xlim([0 18]) ; % ylim([0 0.04]) ;
title('Single')
ylabel('MT_{sat}); xlabel('B_{1rms}');
ax = gca; ax.FontSize = 20; 
colorbar off
hold off                     
  set(gcf,'position',[10,400,1200,400])                       
     sgtitle(str)             
     
 
saveas(gcf,strcat(SavDir,'initial_best_parameters_fit2Optimal_2.png'))  
     







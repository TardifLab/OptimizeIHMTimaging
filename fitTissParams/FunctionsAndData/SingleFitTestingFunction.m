function idx = SingleFitTestingFunction(baseDir, saveImgDir, TParams, idx)


%% Load data
ref_kspace_dir = strcat(baseDir,'kspaceWeighting/Atlas_reference/');
gm_m = load( strcat( ref_kspace_dir,'GM_seg_MNI_152_image.mat')); % keep these here so you know what files you are setting directory to
fft_gm_m = load( strcat( ref_kspace_dir,'GM_seg_MNI_152_kspace.mat'));

gm_m = gm_m.gm_m;
fft_gm_m = fft_gm_m.fft_gm_m;

exportMat = load( [baseDir,'fitTissParams/MTsat_vals2fit.mat'] );   
exportMat2 = load( [baseDir,'fitTissParams/MTsat_vals2fit_Mar16.mat'] );   

exportMat = exportMat.exportMat;
exportMat2 = exportMat2.exportMat2;


t1 = TParams(1);
t2 = TParams(2);
t3 = TParams(3);
t4 = TParams(4);
t5 = TParams(5);
t6 = TParams(6);

Single_sig1_f = zeros(  Params.TurboFactor, length(B1rms));
Single_sig2_f = zeros(  Params2.TurboFactor, length(B1rms));
Single_sig3_f = zeros(  Params3.TurboFactor, length(B1rms));
Dual_sig1_f = zeros( Params.TurboFactor,  length(B1rms));
Dual_sig2_f = zeros( Params2.TurboFactor, length(B1rms));
Dual_sig3_f = zeros( Params3.TurboFactor, length(B1rms));


parfor i = 1:10
     [Ssig1,~, ~]  = BlochSimFlashSequence_v2(Params,'freqPattern', 'single', 'b1', B1rms(i),...
         'R', t1, 'T2a', t2, 'T1D',t3, 'T2b', t4, 'M0b', t5, 'R1b', t6 );           

     [Ssig2,~, ~]  = BlochSimFlashSequence_v2(Params2,'freqPattern', 'single', 'b1', B1rms(i),...
         'R', t1, 'T2a', t2, 'T1D',t3, 'T2b', t4, 'M0b', t5, 'R1b', t6 );

     [Ssig3,~, ~]  = BlochSimFlashSequence_v2(Params3,'freqPattern', 'single', 'b1', B1rms(i),...
         'R', t1, 'T2a', t2, 'T1D',t3, 'T2b', t4, 'M0b', t5, 'R1b', t6 );

     [Dsig1,~, ~]  = BlochSimFlashSequence_v2(Params,'freqPattern', 'dualAlternate', 'b1', B1rms(i),...
         'R', t1, 'T2a', t2, 'T1D',t3, 'T2b', t4, 'M0b', t5, 'R1b', t6 );

     [Dsig2,~, ~]  = BlochSimFlashSequence_v2(Params2,'freqPattern', 'dualAlternate', 'b1', B1rms(i),...
         'R', t1, 'T2a', t2, 'T1D',t3, 'T2b', t4, 'M0b', t5, 'R1b', t6 );

     [Dsig3,~, ~]  = BlochSimFlashSequence_v2(Params3,'freqPattern', 'dualAlternate', 'b1', B1rms(i),...
         'R', t1, 'T2a', t2, 'T1D',t3, 'T2b', t4, 'M0b', t5, 'R1b', t6 ); 

     % store values
     Single_sig1_f(:,i)   = Ssig1;
     Single_sig2_f(:,i)   = Ssig2;
     Single_sig3_f(:,i)   = Ssig3;

     Dual_sig1_f(:,i)   = Dsig1;
     Dual_sig2_f(:,i)   = Dsig2;
     Dual_sig3_f(:,i)   = Dsig3;
end


Single_sig1_gm_f = zeros(length(B1rms),1);
Single_sig2_gm_f = zeros(length(B1rms),1);
Single_sig3_gm_f = zeros(length(B1rms),1);
Dual_sig1_gm_f = zeros(length(B1rms),1);
Dual_sig2_gm_f = zeros(length(B1rms),1);
Dual_sig3_gm_f = zeros(length(B1rms),1);

% calculate values based on full k-space data
for j = 1:length(B1rms)
    Single_sig1_gm_f(j) = CR_generate_BSF_scaling_v1( squeeze(Single_sig1_f(:,j)), Params, outputSamplingTable, gm_m, fft_gm_m) ;  
    Single_sig2_gm_f(j) = CR_generate_BSF_scaling_v1(squeeze(Single_sig2_f(:,j)), Params2, outputSamplingTable2, gm_m, fft_gm_m) ;
    Single_sig3_gm_f(j) = CR_generate_BSF_scaling_v1(squeeze(Single_sig3_f(:,j)), Params3, outputSamplingTable3, gm_m, fft_gm_m) ;
    Dual_sig1_gm_f(j) = CR_generate_BSF_scaling_v1(squeeze(Dual_sig1_f(:,j)), Params, outputSamplingTable, gm_m, fft_gm_m) ;
    Dual_sig2_gm_f(j) = CR_generate_BSF_scaling_v1(squeeze(Dual_sig2_f(:,j)), Params2, outputSamplingTable2, gm_m, fft_gm_m);
    Dual_sig3_gm_f(j) = CR_generate_BSF_scaling_v1(squeeze(Dual_sig3_f(:,j)), Params3, outputSamplingTable3, gm_m, fft_gm_m);   
end


% calculate MTsat
T1obs = ones(size(Dual_sig1_gm_f)) .* 1.4.*1000;
M0_app_v = ones(size(Dual_sig1_gm_f)) ;


MTsat_sim_Single1_f = calcMTsatThruLookupTablewithDummyV3( Single_sig1_gm_f, [], T1obs, [], M0_app_v, Params.echoSpacing * 1000,  Params.numExcitation,  Params.TR * 1000,  Params.flipAngle,   Params.DummyEcho);
MTsat_sim_Single2_f = calcMTsatThruLookupTablewithDummyV3( Single_sig2_gm_f, [], T1obs, [], M0_app_v, Params2.echoSpacing * 1000, Params2.numExcitation, Params2.TR * 1000, Params2.flipAngle, Params2.DummyEcho);
MTsat_sim_Single3_f = calcMTsatThruLookupTablewithDummyV3( Single_sig3_gm_f, [], T1obs, [], M0_app_v, Params3.echoSpacing * 1000, Params3.numExcitation, Params3.TR * 1000, Params3.flipAngle, Params3.DummyEcho);
MTsat_sim_Dual1_f   = calcMTsatThruLookupTablewithDummyV3( Dual_sig1_gm_f,   [], T1obs, [], M0_app_v, Params.echoSpacing * 1000,  Params.numExcitation,  Params.TR * 1000,  Params.flipAngle,   Params.DummyEcho);
MTsat_sim_Dual2_f   = calcMTsatThruLookupTablewithDummyV3( Dual_sig2_gm_f,   [], T1obs, [], M0_app_v, Params2.echoSpacing * 1000, Params2.numExcitation, Params2.TR * 1000, Params2.flipAngle, Params2.DummyEcho);
MTsat_sim_Dual3_f   = calcMTsatThruLookupTablewithDummyV3( Dual_sig3_gm_f,   [], T1obs, [], M0_app_v, Params3.echoSpacing * 1000, Params3.numExcitation, Params3.TR * 1000, Params3.flipAngle, Params3.DummyEcho);

% fit the simulations with a polynomial
Single_c1_f = polyfit(B1rms, MTsat_sim_Single1_f, fit_degree);
Single_c2_f = polyfit(B1rms, MTsat_sim_Single2_f, fit_degree);
Single_c3_f = polyfit(B1rms, MTsat_sim_Single3_f, fit_degree);
Dual_c1_f   = polyfit(B1rms, MTsat_sim_Dual1_f, fit_degree);
Dual_c2_f   = polyfit(B1rms, MTsat_sim_Dual2_f, fit_degree);
Dual_c3_f   = polyfit(B1rms, MTsat_sim_Dual3_f, fit_degree);

% combine both scatterplots


[~, sortidx_b1] = sort(  exportMat(10,:) , 'ascend');

% Sort the whole matrix:
mat_sort = exportMat( :, sortidx_b1);

% Smooth matrix
mat_ss1 = smoothdata( mat_sort,2, 'movmedian', 10);

% Make sure you have B1rms values! Multiply B1 map by the B1rms of the sequence 
b1_1_1 =11.4* mat_ss1(10,:) ; 
b1_2_1 = 13.3 * mat_ss1(10,:) ; 
b1_3_1 = 11.6*mat_ss1(10,:) ; 


% stack matrices:
% Sort B1 , then sat values, then smooth them
[~, sortidx_b1] = sort(  exportMat2(10,:) , 'ascend');

% Sort the whole matrix:
mat_sort = exportMat2( :, sortidx_b1);

% Smooth matrix
mat_ss2 = smoothdata( mat_sort,2, 'movmedian', 10);

% Make sure you have B1rms values! Multiply B1 map by the B1rms of the sequence 
b1_1_2 =11.4* mat_ss2(10,:) ; 
b1_2_2 = 13.3 * mat_ss2(10,:) ; 
b1_3_2 = 11.6*mat_ss2(10,:) ; 


% concatenate:
b1_1 = [b1_1_1, b1_1_2];
b1_2 = [b1_2_1, b1_2_2];
b1_3 = [b1_3_1, b1_3_2];

mat_ss = [mat_ss1, mat_ss2];



%% Now plot the results
x1_line = linspace(0,18,100);
y1 = polyval( Dual_c1_f, x1_line);
y2 = polyval( Dual_c2_f, x1_line);
y3 = polyval( Dual_c3_f, x1_line);

str = ['R = ',num2str(TParams(1)),', T2a = ',num2str(TParams(2)),...
    ', T1D = ',num2str(TParams(3)), ', T2b =',num2str(TParams(4)),...
    ', M0b = ',num2str(TParams(5)),', R1b = ',num2str(TParams(6))];

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


y1 = polyval( Single_c1_f, x1_line);
y2 = polyval( Single_c2_f, x1_line);
y3 = polyval( Single_c3_f, x1_line);

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
        

saveas(gcf,strcat(saveImgDir,'Final_best_parameters_fit2Optimal_',num2str(idx),'.png'))  
idx = idx+1;    
     

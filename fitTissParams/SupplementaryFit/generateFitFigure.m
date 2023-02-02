function generateFitFigure(Params, Params2, Params3,mat_ss,...
        M0B, R,  T2B,  T1D, T2A, R1B,Savefn,...
        gm_m, fft_gm_m,outputSamplingTable, outputSamplingTable2, outputSamplingTable3 )

B1rms = linspace(0,1.3,14);

Ssig1 = zeros(Params.TurboFactor, length(B1rms) );
Ssig2 = zeros(Params2.TurboFactor, length(B1rms) );
Ssig3 = zeros(Params3.TurboFactor, length(B1rms) );
Dsig1 = zeros(Params.TurboFactor, length(B1rms) );
Dsig2 = zeros(Params2.TurboFactor, length(B1rms) );
Dsig3 = zeros(Params3.TurboFactor, length(B1rms) );

% simulate for each sequence, length(B1rms) relative b1 values, get mag over readout
parfor i = 1:length(B1rms)
     [Ssig1(:,i),~, ~]  = BlochSimFlashSequence_v2(Params,'freqPattern', 'single', 'b1', Params.b1*B1rms(i),...
         'R', R, 'T2a', T2A, 'T1D',T1D, 'T2b', T2B, 'M0b', M0B, 'R1b', R1B ); 
                    
     [Ssig2(:,i),~, ~]  = BlochSimFlashSequence_v2(Params2,'freqPattern', 'single', 'b1', Params2.b1*B1rms(i),...
         'R', R, 'T2a', T2A, 'T1D',T1D, 'T2b', T2B, 'M0b', M0B, 'R1b', R1B );
     
     [Ssig3(:,i),~, ~]  = BlochSimFlashSequence_v2(Params3,'freqPattern', 'single', 'b1', Params3.b1*B1rms(i),...
         'R', R, 'T2a', T2A, 'T1D',T1D, 'T2b', T2B, 'M0b', M0B, 'R1b', R1B );
     
     [Dsig1(:,i),~, ~]  = BlochSimFlashSequence_v2(Params,'freqPattern', 'dualAlternate', 'b1', Params.b1* B1rms(i),...
         'R', R, 'T2a', T2A, 'T1D',T1D, 'T2b', T2B, 'M0b', M0B, 'R1b', R1B );
          
     [Dsig2(:,i),~, ~]  = BlochSimFlashSequence_v2(Params2,'freqPattern', 'dualAlternate', 'b1', Params2.b1*B1rms(i),...
         'R', R, 'T2a', T2A, 'T1D',T1D, 'T2b', T2B, 'M0b', M0B, 'R1b', R1B );

     [Dsig3(:,i),~, ~]  = BlochSimFlashSequence_v2(Params3,'freqPattern', 'dualAlternate', 'b1', Params3.b1* B1rms(i),...
         'R', R, 'T2a', T2A, 'T1D',T1D, 'T2b', T2B, 'M0b', M0B, 'R1b', R1B ); 
end

%% convert readout to image intensity:
sig = zeros( length(B1rms), 6);

for i = 1:length(B1rms) % for each b1
    sig(i,1) = CR_generate_BSF_scaling_v1( squeeze(Ssig1(:,i)), Params,  outputSamplingTable,  gm_m, fft_gm_m) ;  
    sig(i,2) = CR_generate_BSF_scaling_v1( squeeze(Ssig2(:,i)), Params2, outputSamplingTable2, gm_m, fft_gm_m) ;
    sig(i,3) = CR_generate_BSF_scaling_v1( squeeze(Ssig3(:,i)), Params3, outputSamplingTable3, gm_m, fft_gm_m) ;
    sig(i,4) = CR_generate_BSF_scaling_v1( squeeze(Dsig1(:,i)), Params,  outputSamplingTable,  gm_m, fft_gm_m) ;
    sig(i,5) = CR_generate_BSF_scaling_v1( squeeze(Dsig2(:,i)), Params2, outputSamplingTable2, gm_m, fft_gm_m) ;
    sig(i,6) = CR_generate_BSF_scaling_v1( squeeze(Dsig3(:,i)), Params3, outputSamplingTable3, gm_m, fft_gm_m) ;  
end

%% Use intensity to calculate MTsat
sat = zeros( length(B1rms), 6);
T1obs = ones(length(B1rms),1) .* 1.4.*1000;
M0_app = ones(length(B1rms),1);

sat(:,1) = calcMTsatThruLookupTablewithDummyV3( sig(:,1), [], T1obs, [], M0_app, Params.echoSpacing * 1000,  Params.numExcitation,  Params.TR * 1000,  Params.flipAngle,   Params.DummyEcho);
sat(:,2) = calcMTsatThruLookupTablewithDummyV3( sig(:,2), [], T1obs, [], M0_app, Params2.echoSpacing * 1000, Params2.numExcitation, Params2.TR * 1000, Params2.flipAngle, Params2.DummyEcho);
sat(:,3) = calcMTsatThruLookupTablewithDummyV3( sig(:,3), [], T1obs, [], M0_app, Params3.echoSpacing * 1000, Params3.numExcitation, Params3.TR * 1000, Params3.flipAngle, Params3.DummyEcho);
sat(:,4) = calcMTsatThruLookupTablewithDummyV3( sig(:,4), [], T1obs, [], M0_app, Params.echoSpacing * 1000,  Params.numExcitation,  Params.TR * 1000,  Params.flipAngle,   Params.DummyEcho);
sat(:,5) = calcMTsatThruLookupTablewithDummyV3( sig(:,5), [], T1obs, [], M0_app, Params2.echoSpacing * 1000, Params2.numExcitation, Params2.TR * 1000, Params2.flipAngle, Params2.DummyEcho);
sat(:,6) = calcMTsatThruLookupTablewithDummyV3( sig(:,6), [], T1obs, [], M0_app, Params3.echoSpacing * 1000, Params3.numExcitation, Params3.TR * 1000, Params3.flipAngle, Params3.DummyEcho);

%% Fit the SIMULATION data

fit_degree = 7;

Single_c1 = polyfit(B1rms, sat(:,1), fit_degree);
Single_c2 = polyfit(B1rms, sat(:,2), fit_degree);
Single_c3 = polyfit(B1rms, sat(:,3), fit_degree);
Dual_c1   = polyfit(B1rms, sat(:,4), fit_degree);
Dual_c2   = polyfit(B1rms, sat(:,5), fit_degree);
Dual_c3   = polyfit(B1rms, sat(:,6), fit_degree);


str = ['M0B = ',num2str(M0B),', R = ',num2str(R),...
    ', T2B = ',num2str(T2B), ', T1D =',num2str(T1D),...
    ', T2A = ',num2str(T2A),', R1B = ',num2str(R1B)];
            

x1_line = linspace(0,1.3,50);
y1 = polyval( Dual_c1, x1_line);
y2 = polyval( Dual_c2, x1_line);
y3 = polyval( Dual_c3, x1_line);
b1 = mat_ss(10,:)';

figure;
subplot(1,2,1)
heatscatter(b1, mat_ss(1,:)' ); 
hold on
heatscatter(b1, mat_ss(2,:)' ); 
heatscatter(b1, mat_ss(3,:)' ); 
hold on
plot(x1_line,y1,'LineWidth',3); plot(x1_line,y2,'LineWidth',3); plot(x1_line,y3,'LineWidth',3)
xlim([0 1.3]) ; % ylim([0 0.04]) ;
title('Dual');
colorbar off
ax = gca; ax.FontSize = 20; 
hold off


y1 = polyval( Single_c1, x1_line);
y2 = polyval( Single_c2, x1_line);
y3 = polyval( Single_c3, x1_line);

subplot(1,2,2)
heatscatter(b1, mat_ss(4,:)' ); 
hold on
heatscatter(b1, mat_ss(5,:)' ); 
heatscatter(b1, mat_ss(6,:)' ); 
hold on
plot(x1_line,y1,'LineWidth',3); plot(x1_line,y2,'LineWidth',3); plot(x1_line,y3,'LineWidth',3)
xlim([0 1.3]) ; % ylim([0 0.04]) ;
title('Single')
ax = gca; ax.FontSize = 20; 
colorbar off
hold off                     
  set(gcf,'position',[10,400,1200,400])                       
     sgtitle(str)             
     
     
   saveas(gcf,Savefn)  
     












function sat = simSequenceForFit(Params, Params2, Params3,...
        M0B, R,  T2B,  T1D, T2A, R1B, B1rms,...
        gm_m, fft_gm_m,outputSamplingTable, outputSamplingTable2, outputSamplingTable3 )

% This simulates the signal, and outputs a matrix where each row is a
% separate B1 relative value, with six columns for the 6 different
% sequences. outputs are the MTsat values
Dsig1 = zeros(Params.TurboFactor, 4 );
Dsig2 = zeros(Params2.TurboFactor, 4 );
Dsig3 = zeros(Params3.TurboFactor, 4 );
Ssig1 = zeros(Params.TurboFactor, 4 );
Ssig2 = zeros(Params2.TurboFactor, 4 );
Ssig3 = zeros(Params3.TurboFactor, 4 );

% simulate for each sequence, 4 relative b1 values, get mag over readout
for i = 1:4
     [Dsig1(:,i),~, ~]  = BlochSimFlashSequence_v2(Params,'freqPattern', 'dualAlternate', 'b1', Params.b1*B1rms(i),...
         'R', R, 'T2a', T2A, 'T1D',T1D, 'T2b', T2B, 'M0b', M0B, 'R1b', R1B ); 
                    
     [Dsig2(:,i),~, ~]  = BlochSimFlashSequence_v2(Params2,'freqPattern', 'dualAlternate', 'b1', Params2.b1*B1rms(i),...
         'R', R, 'T2a', T2A, 'T1D',T1D, 'T2b', T2B, 'M0b', M0B, 'R1b', R1B );
     
     [Dsig3(:,i),~, ~]  = BlochSimFlashSequence_v2(Params3,'freqPattern', 'dualAlternate', 'b1', Params3.b1*B1rms(i),...
         'R', R, 'T2a', T2A, 'T1D',T1D, 'T2b', T2B, 'M0b', M0B, 'R1b', R1B );
     
     [Ssig1(:,i),~, ~]  = BlochSimFlashSequence_v2(Params,'freqPattern', 'single', 'b1', Params.b1* B1rms(i),...
         'R', R, 'T2a', T2A, 'T1D',T1D, 'T2b', T2B, 'M0b', M0B, 'R1b', R1B );
          
     [Ssig2(:,i),~, ~]  = BlochSimFlashSequence_v2(Params2,'freqPattern', 'single', 'b1', Params2.b1*B1rms(i),...
         'R', R, 'T2a', T2A, 'T1D',T1D, 'T2b', T2B, 'M0b', M0B, 'R1b', R1B );

     [Ssig3(:,i),~, ~]  = BlochSimFlashSequence_v2(Params3,'freqPattern', 'single', 'b1', Params3.b1* B1rms(i),...
         'R', R, 'T2a', T2A, 'T1D',T1D, 'T2b', T2B, 'M0b', M0B, 'R1b', R1B ); 
end

%% convert readout to image intensity:
sig = zeros( 4, 6);

for i = 1:4 % for each b1
    sig(i,1) = CR_generate_BSF_scaling_v1( squeeze(Dsig1(:,i)), Params,  outputSamplingTable,  gm_m, fft_gm_m) ;  
    sig(i,2) = CR_generate_BSF_scaling_v1( squeeze(Dsig2(:,i)), Params2, outputSamplingTable2, gm_m, fft_gm_m) ;
    sig(i,3) = CR_generate_BSF_scaling_v1( squeeze(Dsig3(:,i)), Params3, outputSamplingTable3, gm_m, fft_gm_m) ;
    sig(i,4) = CR_generate_BSF_scaling_v1( squeeze(Ssig1(:,i)), Params,  outputSamplingTable,  gm_m, fft_gm_m) ;
    sig(i,5) = CR_generate_BSF_scaling_v1( squeeze(Ssig2(:,i)), Params2, outputSamplingTable2, gm_m, fft_gm_m) ;
    sig(i,6) = CR_generate_BSF_scaling_v1( squeeze(Ssig3(:,i)), Params3, outputSamplingTable3, gm_m, fft_gm_m) ;  
end

%% Use intensity to calculate MTsat
sat = zeros( 4, 6);
T1obs = ones(4,1) .* 1.4.*1000;
M0_app = ones(4,1);

sat(:,1) = calcMTsatThruLookupTablewithDummyV3( sig(:,1), [], T1obs, [], M0_app, Params.echoSpacing * 1000,  Params.numExcitation,  Params.TR * 1000,  Params.flipAngle,   Params.DummyEcho);
sat(:,2) = calcMTsatThruLookupTablewithDummyV3( sig(:,2), [], T1obs, [], M0_app, Params2.echoSpacing * 1000, Params2.numExcitation, Params2.TR * 1000, Params2.flipAngle, Params2.DummyEcho);
sat(:,3) = calcMTsatThruLookupTablewithDummyV3( sig(:,3), [], T1obs, [], M0_app, Params3.echoSpacing * 1000, Params3.numExcitation, Params3.TR * 1000, Params3.flipAngle, Params3.DummyEcho);
sat(:,4) = calcMTsatThruLookupTablewithDummyV3( sig(:,4), [], T1obs, [], M0_app, Params.echoSpacing * 1000,  Params.numExcitation,  Params.TR * 1000,  Params.flipAngle,   Params.DummyEcho);
sat(:,5) = calcMTsatThruLookupTablewithDummyV3( sig(:,5), [], T1obs, [], M0_app, Params2.echoSpacing * 1000, Params2.numExcitation, Params2.TR * 1000, Params2.flipAngle, Params2.DummyEcho);
sat(:,6) = calcMTsatThruLookupTablewithDummyV3( sig(:,6), [], T1obs, [], M0_app, Params3.echoSpacing * 1000, Params3.numExcitation, Params3.TR * 1000, Params3.flipAngle, Params3.DummyEcho);





















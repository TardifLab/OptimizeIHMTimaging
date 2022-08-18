function simResults = CR_add_ihMTsat_SNR(simResults, ParameterSet, noiseLvl, T1, MP2RAGE)
% to be used with plotSimResultsIHMT_figures_SNR_v1.m
% Requires MP2RAGE scripts for this: https://github.com/JosePMarques/MP2RAGE-related-scripts

if nargin < 3 || isempty(noiseLvl)
  noiseLvl = 0.0005; 
end
if nargin < 4 || isempty(T1)
  T1 = 1.4;
end
if nargin < 5 || isempty(MP2RAGE)
    MP2RAGE.B0=3;           % in Tesla
    MP2RAGE.TR=5;           % MP2RAGE TR in seconds 
    MP2RAGE.TRFLASH=6.4e-3; % TR of the GRE readout
    MP2RAGE.TIs=[940e-3 2830e-3];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
    MP2RAGE.NZslices=176;% C.R. just need the turbofactor here.
    MP2RAGE.FlipDegrees=[4 5];% Flip angle of the two readouts in degrees
end


% Generate noise vectors (30,000 x 4)
additive_noise = normrnd( 1,noiseLvl ,30000,5) - 1;
%std(additive_noise) % Below 30k, you don't get the same value

% can check visually the result:
%histogram(additive_noise(:,1))

% Use the T1 and M0 == 1 to get simulated signal values for MP2RAGE
MP2RAGE_Signal = 1*MPRAGEfunc(2, MP2RAGE.TR, MP2RAGE.TIs, MP2RAGE.NZslices, MP2RAGE.TRFLASH, MP2RAGE.FlipDegrees, 'normal',T1);

% Make 2 vectors, 1 for each inversion time, and add the random noise
inv1.img = additive_noise(:,1) + MP2RAGE_Signal(1);
inv2.img = additive_noise(:,2) + MP2RAGE_Signal(2);

% Generate Uni Image
Uni.img=real(inv1.img).*conj(inv2.img)./(abs(inv1.img).^2+abs(inv2.img).^2); % taken from line 68 of MP2RAGE_lookuptable.m
B1.img = ones(size(inv2.img));

% uncenter and multiply up
Uni.img = double(Uni.img+0.5)*4095;

% Calculate noisy vectors for M0 and T1
[ T1_n, ~, M0_n] = CR_T1B1correctpackageTFL_withM0( B1, Uni, inv2, MP2RAGE, [], 1);

T1_n = T1_n.img; % want T1 in ms for lookup table
M0_n = M0_n.img;

% View results if desired
% histogram(M0_n)
% histogram(T1_n)

%% For MTsat
echospacing = 7.66; % in ms
ihMTsat_SNR = zeros(size(simResults(:,1)));

noise1 = additive_noise(:,3);
noise2 = additive_noise(:,4);
noise3 = additive_noise(:,5);

%% Loop this bit for each row:
parfor i = 1:size(simResults,1)

    % Generate noisy signals
    Single_n = simResults(i,1) + noise1;
    Single_n2 = simResults(i,1) + noise2;
    Dual_n   = simResults(i,2) + noise3;
    
    if ParameterSet(i,6) == 1
        DummyEcho = 0;
    elseif ParameterSet(i,6) < 4
        DummyEcho = 1;
    else 
        DummyEcho = 2;
    end

    
    MTsat_s_n = calcMTsatThruLookupTablewithDummyV3( Single_n, [],T1_n, ...
        [], M0_n, echospacing, ParameterSet(i,6) + DummyEcho,...
        ParameterSet(i,3)*1000, ParameterSet(i,2), DummyEcho);


    MTsat_sn2 = calcMTsatThruLookupTablewithDummyV3( Single_n2, [],T1_n, ...
        [], M0_n, echospacing, ParameterSet(i,6) + DummyEcho,...
        ParameterSet(i,3)*1000, ParameterSet(i,2), DummyEcho);


    MTsat_d_n = calcMTsatThruLookupTablewithDummyV3( Dual_n, [],T1_n,...
        [], M0_n, echospacing, ParameterSet(i,6) + DummyEcho,...
        ParameterSet(i,3)*1000, ParameterSet(i,2), DummyEcho);
    

    % Noisy ihMTsat
    ihMTsat_n = MTsat_d_n - (MTsat_s_n+MTsat_sn2)./2;

    ihMT_noise = std(ihMTsat_n);
    ihMTsat_SNR(i) = mean(ihMTsat_n)./ ihMT_noise;
end



ihMT_SNR_eff = ihMTsat_SNR ./simResults(:,10); % efficiency, where acquisition time is in column 10.

simResults(:,13) = ihMTsat_SNR;
simResults(:,14) = ihMT_SNR_eff;

%     histogram(MTsat_s_n)
%     histogram(MTsat_d_n)
%     histogram(ihMTsat_n)






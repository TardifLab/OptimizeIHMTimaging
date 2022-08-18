%% This code is meant to solve the signal equations presented in Deichmann et al., 2000 for optimizing an MP-RAGE sequence
% The goal of this project is a multi-segment gre readout for ihMT imaging.
function Sig= MTrage_sig_eqn_withDummy(echospacing, flip, T1, TD, numReadExcitation, M0, MT_drop, B1field, MT_b1_corr, dummyEchoes)

% NOTE: here numReadExcitation == turbofactor - dummyEchoes

%% Test parameters
% echospacing= Params.echospacing; % The echo spacing of the GRE readout
% numExcitation = Params.numExcitation; % in my implementation, only odd numbers work. 
% TD = Params.TD; % dead time in sequence, for SAR 
% flip = Params.flipAngle; % flip angle
% MT_drop = 0; % At steady state, how much does the signal drop? 
% B1field = 1;
% MT_b1_corr = 1;
% T1 = 1;
% M0 = 1;
% dummyEchoes = 2  = number of dummy echoes at start of each echo train

%% Equations
flip_a = (flip*B1field) * pi / 180; % correct for B1 and convert to radians

%% Following readout you magnetization (M2) = A1 + M1 * B , derivation at the bottom of function. 

x = cos(flip_a) ;
y = exp(-echospacing/T1);

B1 = (x*y)^numReadExcitation;
% A1 = 0;
% for i = 1:numExcitation
%    A1 = A1 + M0*(x*y)^(i-1);
%    A1 = A1 - M0*(x*y)^(i)/x ;
% end

%% redo based on result from Munsch et al 2021 ihMT paper:
A1 = M0*(1-y)* ( (1- x^numReadExcitation * exp(-numReadExcitation*echospacing/T1)) / (1- x*y) );

%% You then have some time TD for T1 relaxation, M3 = A2 + B2*M2
A2 = M0 * (1 - exp(-TD / T1));
B2 =  exp( -TD  / T1    );

%% Following TD, your MTsat pulses knock down your longitudinal relaxation by factor MTdrop M4 = A3 + M3 * B3
A3 = 0; % we will assume it is instanteous, and the relaxation time is in the one above.
B3 = (1-(MT_drop*MT_b1_corr));

%% You might have some dummy echoes played out before long echo train. M5 = A4 + M4*B4

A4 = M0*(1-y)* ( (1- x^dummyEchoes * exp(-dummyEchoes*echospacing/T1)) / (1- x*y) );
B4 = (x*y)^dummyEchoes;


%% Since M5 = M1 and we want to solve for M1
M = (A4 + A3*B4 + A2*B3*B4 + A1*B2*B3*B4) / (1-B1*B2*B3*B4); % + A1*B2*B3

% With central encoding we readout M with sin(flip)
Sig = sin(flip_a) * M;




% % DERIVE M1 = A/(1-B) 
% M5 = A4 + B4 * M4
% M5 = A4 + B4 * (A3 + B3 * M3)
% M5 = A4 + B4 * (A3 + B3 * (A2 + B2 * (A1 + B1 * M1))
% M5 = A4 + B4 * (A3 + B3 * (A2 + A1*B2 +B1*B2*M1))
% M5 = A4 + B4 * (A3 + A2*B3 + A1*B2*B3 +B1*B2*B3*M1)
% M5 = A4 + A3*B4 + A2*B3*B4 + A1*B2*B3*B4 +B1*B2*B3*B4*M1 % DIVIDE BY M1 (==M5)
% 1 = (A4 + A3*B4 + A2*B3*B4 + A1*B2*B3*B4 +B1*B2*B3*B4*M1)/M1
% M1 - B1*B2*B3*B4*M1 = A4 + A3*B4 + A2*B3*B4 + A1*B2*B3*B4
% 
% Set A = A4 + A3*B4 + A2*B3*B4 + A1*B2*B3*B4
% Set B = B3*B2*B1
% 
% M1 - B*M1 = A
% M1(1-B) = A
% M1 = A/(1-B)
















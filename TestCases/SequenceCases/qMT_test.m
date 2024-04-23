% Test to see how the simulation code stacks against the theoretical value
% from a basic flash sequence for VFA T1 and M0 mapping


% Use the optimizations from: Levesque, I.R., Sled, J.G., Pike, G.B., 2011. 
% Iterative optimization method for design of quantitative magnetization transfer 
% imaging experiments. Magn. Reson. Med. 66, 635–643. https://doi.org/10.1002/mrm.23071

function qMT_resultsTable = qMT_test(savDir)

%savDir =  'C:\Users\crowle1\OneDrive - McGill University\ihMT_work\cortical_ihMT_sim\simCode\sim_wSpoil\RF_grad_diffusion_v4\TestCases\Figures\';
Params.SatPulseShape = 'gausshann';

%% Note to compare to their results, use 1.5T
Params.B0 = 1.5;
Params.MTC = 0; % Magnetization Transfer Contrast
Params.TissueType = 'WM';
Params.echoSpacing = 5/1000;
Params.numExcitation = 1;

Params = DefaultCortexTissueParams(Params);
Params = CalcImagingParams(Params);
Params = CalcVariableImagingParams(Params);
Params.CalcVector = 1;


%% Generate a list of qMT parameters:
pulseDur1 = 10.24/1000;
pulseDur2 = 30.72/1000;

% Copying table 1 for 1.5T
% [ TR, flip angle, alpha(MT), delta (Hz), relSig  (estimated from plot 3a.)
% In that plot, solid line is 25ms TR, dotted = 50ms.
MTParams = [...
    25,  7, 142, 800,   0.84;...
    25,  7, 284, 2010,  0.72; ...
    25,  7, 426, 9372,  0.85; ...
    25,  7, 710, 1478,  0.42; ...
    25,  7, 568, 12679, 0.88; ...
    50, 10, 347, 127,   0.46; ...
    25,  7, 710, 1088,  0.39; ...
    50, 10, 694, 2010,  0.72; ...
    50, 10, 347, 172,   0.56; ...
    25,  7, 710, 12679, 0.83]; %; ... % plot only shows best 10
%     25,  7, 710, 2010, 0   ; ...
%     50, 10, 694, 1478, 0    ; ...
%     50, 10, 347, 94, 0      ; ...
%     50, 10, 694, 2732, 0    ; ...
%     25,  7, 710, 800, 0]    ;

Prot1 = MTParams(:,1) == 25;
Prot2 = MTParams(:,1) == 50;

% convert TR to seconds.
MTParams(:,1) = MTParams(:,1)/1000; 

% % view table plotted:
% figure;
% scatter(MTParams(:,4),MTParams(:,5),20,'r','filled')
% set(gca, 'XScale', 'log'); grid on
% ylim([0 1]); xlim([10 10^5])


%% Sim sequences to get resulting values
% they do prep pulses for 12seconds, I am seeing steady state with these
% pulses before 5 seconds, so I won't modify my code. 

% ref val 25ms
Params.TR = 25/1000;
Params.flipAngle = 7;
Params.MTC = 0; % Magnetization Transfer Contrast

[ref25, ~, ~] = BlochSimFlashSequence_v2(Params, 'GradientSpoiling',1,'ModelSpinDiffusion', 1);
%[ref25, ~, ~] = BlochSimFlashSequence_v2(Params, 'PerfectSpoiling',1,'IncludeDipolar',0);

% ref val 50ms
Params.TR = 50/1000;
Params.flipAngle = 10;

[ref50, ~, ~] = BlochSimFlashSequence_v2(Params, 'GradientSpoiling',1,'ModelSpinDiffusion', 1);
%[ref50, ~, ~] = BlochSimFlashSequence_v2(Params, 'PerfectSpoiling',1,'IncludeDipolar',0);

% set a few more variables then simulate
sim_qMT = zeros(1,10);
sim_qMT_norm = zeros(1,10);
Params.MTC = 1;
Params.pulseGapDur = 0.3/1000;
Params.freqPattern = 'single';
Params.SatPulseShape = 'gausshann';
Params.numSatPulse = 1;
Params = CalcVariableImagingParams(Params);
Params.ReadoutResolution = 1.3e-3;

for i = 1:length(sim_qMT)

    Params.TR = MTParams(i,1);
    Params.flipAngle = MTParams(i,2);
    Params.delta = MTParams(i,4);
    Params.satFlipAngle = MTParams(i,3);
    
    if Prot1(i)
        Params.pulseDur = pulseDur1;
        [sim_qMT(i), ~, ~] = BlochSimFlashSequence_v2(Params, 'GradientSpoiling',1,'ModelSpinDiffusion', 1);
        %[sim_qMT(i), ~, ~] = BlochSimFlashSequence_v2(Params, 'PerfectSpoiling',1,'IncludeDipolar',0);
        sim_qMT_norm(i) = sim_qMT(i)/ref25;

    elseif Prot2(i)
        Params.pulseDur = pulseDur2;
        [sim_qMT(i), ~, ~] = BlochSimFlashSequence_v2(Params, 'GradientSpoiling',1,'ModelSpinDiffusion', 1);
        %[sim_qMT(i), ~, ~] = BlochSimFlashSequence_v2(Params, 'PerfectSpoiling',1,'IncludeDipolar',0);
        sim_qMT_norm(i) = sim_qMT(i)/ref50;
    end
end



% Separate for easier interpretation:
P1_paper = MTParams(Prot1,:);
P2_paper = MTParams(Prot2,:);

P1_sim = sim_qMT_norm(Prot1);
P2_sim = sim_qMT_norm(Prot2);

markerSz = 50;

%% Color by b1

clrsc = MTParams(:,3);

mg = mat2gray(clrsc);
mgint = im2uint8(mg);
rgb = squeeze(ind2rgb(mgint, jet));

P1c = rgb(Prot1,:);
P2c = rgb(Prot2,:);

figure;
scatter(P1_paper(:,4),P1_paper(:,5),markerSz,P1c,'o')
hold on
scatter(P2_paper(:,4),P2_paper(:,5),markerSz,P2c,'d')
set(gca, 'XScale', 'log'); grid on
ylim([0 1]); xlim([10 10^5])
scatter(P1_paper(:,4),P1_sim,markerSz,P1c,'x')
scatter(P2_paper(:,4),P2_sim,markerSz,P2c,'+')
legend('Paper TR=25','Paper TR=50','Sim TR=25','Sim TR=50','Location','southeast')
ax = gca; ax.FontSize = 14; 
set(gcf,'position',[10,400,800,600])
xlabel('Frequency Offset (Hz)')
ylabel('Normalized Signal')

saveas(gcf,strcat(savDir,'qMT_test_matchingPoints.png')) 


qMT_resultsTable = table(MTParams(:,5),sim_qMT_norm', 'VariableNames',{'Levesque et al 2011','Simulated'});


save(strcat(savDir,'qMT_resultsTable.mat'), 'qMT_resultsTable')









































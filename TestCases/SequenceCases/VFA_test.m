% Test to see how the simulation code stacks against the theoretical value
% from a basic flash sequence for VFA T1 and M0 mapping

function VFA_resultsTable = VFA_test(savDir)

%savDir =  'C:\Users\crowle1\OneDrive - McGill University\ihMT_work\cortical_ihMT_sim\simCode\sim_wSpoil\RF_grad_diffusion_v4\TestCases\Figures\';


%%  Reference method
Params.TR = 15/1000;
Params.flipAngle = 4;

Params.numExcitation = 1;
Params.MTC = 0; % Magnetization Transfer Contrast
 
Params.B0 = 3;
Params.MTC = 0; % Magnetization Transfer Contrast
Params.TissueType = 'GM';
Params.echoSpacing = 7.7/1000;

Params = DefaultCortexTissueParams(Params);
Params = CalcImagingParams(Params);
Params = CalcVariableImagingParams(Params);
Params.CalcVector = 1;

loops = 400;
M_t0 = zeros(1,2*loops +1);
M_t0(1) = 1;
time_vect = zeros(length(M_t0),1);
idx = 2;
for i = 1:loops
    % Apply Flip Angle
    M_t0(idx) = M_t0(idx-1) * cos (Params.flipAngle  *pi/180);
    time_vect(idx) = time_vect(idx-1);
    idx = idx +1;

    % Relaxation
    M_t0(idx) = 1*(1-exp(-Params.TR*Params.Raobs)) + M_t0(idx-1)*exp(-Params.TR*Params.Raobs);
    time_vect(idx) = time_vect(idx-1) + Params.TR;
    idx = idx +1;
end
M_t0(idx:end) = [];
time_vect(idx:end) = [];

REFERENCE_VALUE1 = M_t0(idx-1)*sin(Params.flipAngle *pi/180);



[lfa0, M1, t1] = BlochSimFlashSequence_v2(Params, 'GradientSpoiling',1,'ModelSpinDiffusion', 0);
[lfa1, M2, t2] = BlochSimFlashSequence_v2(Params, 'GradientSpoiling',1,'ModelSpinDiffusion', 1);
[lfa2, M3, t3] = BlochSimFlashSequence_v2(Params, 'PerfectSpoiling',1);

%% High flip angle:
Params.flipAngle = 20;


M_t1 = zeros(1,2*loops +1);
M_t1(1) = 1;
time_vect1 = zeros(length(M_t1),1);
idx = 2;
for i = 1:loops
    % Apply Flip Angle
    M_t1(idx) = M_t1(idx-1) * cos (Params.flipAngle  *pi/180);
    time_vect1(idx) = time_vect1(idx-1);
    idx = idx +1;

    % Relaxation
    M_t1(idx) = 1*(1-exp(-Params.TR*Params.Raobs)) + M_t1(idx-1)*exp(-Params.TR*Params.Raobs);
    time_vect1(idx) = time_vect1(idx-1) + Params.TR;
    idx = idx +1;
end
M_t1(idx:end) = [];
time_vect1(idx:end) = [];

REFERENCE_VALUE2 = M_t1(idx-1)*sin(Params.flipAngle *pi/180);

[hfa0, M1h, t1h] = BlochSimFlashSequence_v2(Params, 'GradientSpoiling',1,'ModelSpinDiffusion', 0);
[hfa1, M2h, t2h] = BlochSimFlashSequence_v2(Params, 'GradientSpoiling',1,'ModelSpinDiffusion', 1);
[hfa2, M3h, t3h] = BlochSimFlashSequence_v2(Params, 'PerfectSpoiling',1);


%% Calculate metrics
a1 = 4 * pi/180;
a2 = 20 * pi / 180;
TR = Params.TR;

T1_ref =1./( 0.5 .* (REFERENCE_VALUE2.*a2./ TR - REFERENCE_VALUE1.*a1./TR) ./ (REFERENCE_VALUE1./(a1) - REFERENCE_VALUE2./(a2))) *1000;
Aapp_ref = REFERENCE_VALUE1 .* REFERENCE_VALUE2 .* (TR .* a2./a1 - TR.* a1./a2) ./ (REFERENCE_VALUE2.* TR .*a2 - REFERENCE_VALUE1.* TR .*a1);

T1_0 = 1./(0.5 .* (hfa0.*a2./ TR - lfa0.*a1./TR) ./ (lfa0./(a1) - hfa0./(a2))) *1000;
Aapp_0 = lfa0 .* hfa0 .* (TR .* a2./a1 - TR.* a1./a2) ./ (hfa0.* TR .*a2 - lfa0.* TR .*a1);

T1_1 = 1./(0.5 .* (hfa1.*a2./ TR - lfa1.*a1./TR) ./ (lfa1./(a1) - hfa1./(a2))) *1000;
Aapp_1 = lfa1 .* hfa1 .* (TR .* a2./a1 - TR.* a1./a2) ./ (hfa1.* TR .*a2 - lfa1.* TR .*a1);

T1_2 =1./( 0.5 .* (hfa2.*a2./ TR - lfa2.*a1./TR) ./ (lfa2./(a1) - hfa2./(a2))) *1000;
Aapp_2 = lfa2 .* hfa2 .* (TR .* a2./a1 - TR.* a1./a2) ./ (hfa2.* TR .*a2 - lfa2.* TR .*a1);



%% Generate Export Plots:
lg1 = strcat('Z-sim = ',num2str(REFERENCE_VALUE1));
lg2 = strcat('XY no diffusion = ',num2str(lfa0));
lg3 = strcat('XY with diffusion = ',num2str(lfa1));
lg4 = strcat('XY perfect spoil = ',num2str(lfa2));
title_s = strcat('VFA test, flip angle=',num2str( 4 ),', TR=',num2str(Params.TR) );

figure; plot(time_vect, M_t0)
hold on
plot(t1,M1(3,:))
plot(t2,M2(3,:))
plot(t3,M3(3,:))
legend( lg1,lg2,lg3,lg4)
xlabel('Time (s)')
ylabel('M_{za}')
ax = gca; ax.FontSize = 20; 
set(gcf,'position',[10,400,800,600])
title(title_s,'FontSize',20);
xlim([0 5])
saveas(gcf,strcat(savDir,'VFA_test_lfa_mag.png')) 

%%
lg1 = strcat('Z-sim = ',num2str(REFERENCE_VALUE2));
lg2 = strcat('XY no diffusion = ',num2str(hfa0));
lg3 = strcat('XY with diffusion = ',num2str(hfa1));
lg4 = strcat('XY perfect spoil = ',num2str(hfa2));
title_s = strcat('VFA test, flip angle=',num2str( 20 ),', TR=',num2str(Params.TR) );

figure; plot(time_vect1, M_t1)
hold on
plot(t1h,M1h(3,:))
plot(t2h,M2h(3,:))
plot(t3h,M3h(3,:))
legend( lg1,lg2,lg3,lg4)
xlabel('Time (s)')
ylabel('M_{za}')
ax = gca; ax.FontSize = 20; 
set(gcf,'position',[10,400,800,600])
title(title_s,'FontSize',20);
xlim([0 2])
saveas(gcf,strcat(savDir,'VFA_test_hfa_mag.png')) 


%% Build table with results

T1_v = [T1_ref; T1_0; T1_1;T1_2 ];
M0_v = [Aapp_ref; Aapp_0; Aapp_1; Aapp_2];

RowNames = {'Z-sim', 'XY no diffusion','XY with diffusion','XY perfect spoil' };

VFA_resultsTable = table(T1_v,M0_v, 'VariableNames',{'T1(ms)','M0,app'},'RowNames',RowNames);

save(strcat(savDir,'VFA_resultsTable.mat'), 'VFA_resultsTable')






































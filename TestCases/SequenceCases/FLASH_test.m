% Test to see how the simulation code stacks against the theoretical value
% from a basic flash sequence with 2 different parameter sets

function FLASH_test(savDir)

%% Test on a sequence with short TR and higher flip angle
Params.TR = 20/1000;
Params.flipAngle = 20;
Params.numExcitation = 1;
Params.MTC = 0; % Magnetization Transfer Contrast

 
Params.B0 = 3;
Params.MTC = 0; % Magnetization Transfer Contrast
Params.TissueType = 'GM';
Params.echoSpacing = 7.7/1000;

Params = DefaultCortexTissueParams(Params);
Params = CalcImagingParams(Params);
Params = CalcVariableImagingParams(Params);



%%  Reference method
REFERENCE_VALUE1 =  CR_FLASH_solver(Params.flipAngle, Params.TR, 1, Params.M0a, 1./Params.Raobs);

%% Method 1 - Simple Excitation and Relaxation
loops = 500;
M_t2 = zeros(1,2*loops +1);
M_t2(1) = 1;
time_vect1 = zeros(length(M_t2),1);
idx = 2;
for i = 1:loops
    % Apply Flip Angle
    M_t2(idx) = M_t2(idx-1) * cos (Params.flipAngle  *pi/180);
    time_vect1(idx) = time_vect1(idx-1);
    idx = idx +1;

    % Relaxation
    M_t2(idx) = 1*(1-exp(-Params.TR*Params.Raobs)) + M_t2(idx-1)*exp(-Params.TR*Params.Raobs);
    time_vect1(idx) = time_vect1(idx-1) + Params.TR;
    idx = idx +1;
end
M_t2(idx:end) = [];
time_vect1(idx:end) = [];

REFERENCE_VALUE2 = M_t2(idx-1)*sin(Params.flipAngle *pi/180);

%% My simulator:
Params.CalcVector = 1;

tic
[lfa0, M1, t1] = BlochSimFlashSequence_v2(Params, 'GradientSpoiling',1,'ModelSpinDiffusion', 0);

[lfa1, M2, t2] = BlochSimFlashSequence_v2(Params, 'GradientSpoiling',1,'ModelSpinDiffusion', 1);

[lfa2, M3, t3] = BlochSimFlashSequence_v2(Params, 'PerfectSpoiling',1);
toc



% Build legend entries to include resulting signal:
lg1 = strcat('Z-sim = ',num2str(REFERENCE_VALUE2));
lg2 = strcat('XY no diffusion = ',num2str(lfa0));
lg3 = strcat('XY with diffusion = ',num2str(lfa1));
lg4 = strcat('XY perfect spoil = ',num2str(lfa2));
title_s = strcat('FLASH test, flip angle=',num2str(Params.flipAngle),', TR=',num2str(Params.TR) );

figure; plot(time_vect1, M_t2)
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
xlim([0 3])
saveas(gcf,strcat(savDir,'FLASH_test_1.png')) 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Second one with low flip
%% Test on a sequence with short TR and higher flip angle
Params.TR = 20/1000;
Params.flipAngle = 5;
Params.numExcitation = 1;
Params.MTC = 0; % Magnetization Transfer Contrast


%% Method 1 - Simple Excitation and Relaxation
loops = 500;
M_t2 = zeros(1,2*loops +1);
M_t2(1) = 1;
time_vect1 = zeros(length(M_t2),1);
idx = 2;
for i = 1:loops
    % Apply Flip Angle
    M_t2(idx) = M_t2(idx-1) * cos (Params.flipAngle  *pi/180);
    time_vect1(idx) = time_vect1(idx-1);
    idx = idx +1;

    % Relaxation
    M_t2(idx) = 1*(1-exp(-Params.TR*Params.Raobs)) + M_t2(idx-1)*exp(-Params.TR*Params.Raobs);
    time_vect1(idx) = time_vect1(idx-1) + Params.TR;
    idx = idx +1;
end
M_t2(idx:end) = [];
time_vect1(idx:end) = [];

REFERENCE_VALUE2 = M_t2(idx-1)*sin(Params.flipAngle *pi/180);

%% My simulator:
Params.CalcVector = 1;

tic
[lfa0, M1, t1] = BlochSimFlashSequence_v2(Params, 'GradientSpoiling',1,'ModelSpinDiffusion', 0);

[lfa1, M2, t2] = BlochSimFlashSequence_v2(Params, 'GradientSpoiling',1,'ModelSpinDiffusion', 1);

[lfa2, M3, t3] = BlochSimFlashSequence_v2(Params, 'PerfectSpoiling',1);
toc



% Build legend entries to include resulting signal:
lg1 = strcat('Z-sim = ',num2str(REFERENCE_VALUE2));
lg2 = strcat('XY no diffusion = ',num2str(lfa0));
lg3 = strcat('XY with diffusion = ',num2str(lfa1));
lg4 = strcat('XY perfect spoil = ',num2str(lfa2));
title_s = strcat('FLASH test, flip angle=',num2str(Params.flipAngle),', TR=',num2str(Params.TR) );

figure; plot(time_vect1, M_t2)
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
saveas(gcf,strcat(savDir,'FLASH_test_2.png')) 























%% This script is used to look at the point spread function Boosted acquisitions
% Simulations were run to determine best parameters for a range of
% turbofactor. Here we will visualize the impact of the different
% turbofactors on the resulting resolution.

% Written by Christopher Rowley PhD 2022

baseDir = '/path/to/OptimizeIHMTimaging/';

SavDir =  '/path/to/figures/';
SaveImgs = 1; % binary flag for whether to save output figures

turboFactor = [200, 160, 120, 80, 48, 10]; % 

psf_single_values = zeros(length(turboFactor), 21601);
psf_dual_values = zeros(length(turboFactor), 21601);
psf_ihMT_values = zeros(length(turboFactor), 21601);
psf_Ref_values = zeros(length(turboFactor), 21601);

for z = 1:length(turboFactor)
    
     
    %% Step 1 -> simulate the sequence and return the longitudinal magnetization for a full TR
    % This also grabs tissue parameters from DefaultCortexTissueParams()
    Params = CR_getSeqParams( turboFactor(z) ); 

    % Single
    SingleMag = BlochSimFlashSequence_v2(Params);

    % Dual
    Params.freqPattern = 'dualAlternate'; % options: 'single', 'dualAlternate', 'dualContinuous'
    DualMag= BlochSimFlashSequence_v2(Params);

    % Reference for MTR
    Params.b1 = 0; % options: 'single', 'dualAlternate', 'dualContinuous'
    RefMag= BlochSimFlashSequence_v2(Params,'ReferenceScan',1);

    %% Step 2 - fill k-space with table

    Params.NumLines = 216;
    Params.NumPartitions = 192; 
    Params.Slices = 176;
    Params.Grappa = 1;
    Params.ReferenceLines = 32;
    Params.AccelerationFactor = 2;
    Params.Segments = []; 
    %Params.TurboFactor = []; %Params.numExcitation- Params.DummyEcho;
    Params.ellipMask = 1;
    [outputSamplingTable, measuredElem, Params.Segments] = Step1_calculateKspaceSampling_v3 (Params);

%% Can check sampling table with below code
%     figure;
%     imagesc( reshape(outputSamplingTable, Params.NumLines/Params.AccelerationFactor +Params.ReferenceLines, Params.NumPartitions ));
%     caxis([0 200])
%     colormap jet
%     axis image
%     xlabel('Partition')
%     ylabel('Slice')
%     ax = gca;    ax.FontSize = 16;   
%     pause(0.5) % issue with toolbar being in the saved image.
%     if SaveImgs; saveas(gcf,strcat(SavDir,'Kspace_Order_',num2str(Params.TurboFactor),'turbofactor.png')); end

    %% Step 3 -> fill the sampling table with the relative signal values from the simulation, in the order they would be acquired
    kSpaceFill_s = CR_fillKspaceSamplingTable_v2( SingleMag, outputSamplingTable, Params);
    kSpaceFill_d = CR_fillKspaceSamplingTable_v2(   DualMag, outputSamplingTable, Params);
    kSpaceFill_ref = CR_fillKspaceSamplingTable_v2(   RefMag, outputSamplingTable, Params);

    % Add in GRAPPA lines through interpolation
    kSpaceFill_s = CR_interpolateMissingGrappaLines( kSpaceFill_s);
    kSpaceFill_d = CR_interpolateMissingGrappaLines( kSpaceFill_d);
    kSpaceFill_ref = CR_interpolateMissingGrappaLines( kSpaceFill_ref);

    %% Can look at filled sampling table with below code
%     figure; imagesc(kSpaceFill_s); axis image; title('Single Sat'); colorbar; colormap('bone')
%     if SaveImgs; saveas(gcf,strcat(SavDir,'Kspace_Filled_',num2str(Params.TurboFactor),'turbofactor.png')); end
%       pause(0.5)
%     figure; imagesc(kSpaceFill_d); axis image; title('Dual Sat');colorbar; colormap('bone')

%% Step 4 -> Calculate PSF values

    [psf0, interXVec, interYVec] = CR_generate_PSF_upsample( kSpaceFill_s );
    [ psf1, ~,~] = CR_generate_PSF_upsample( kSpaceFill_d );
    [ psf3, ~,~] = CR_generate_PSF_upsample( kSpaceFill_ref );
    psf2 = (psf0 - psf1)./max(psf3,[],'all'); % MTR

    %% Take the center line for an easier display.
    iCenterLine      = floor( size(psf0,1)/2 )+1 ;
    iCenterPartition = floor(  size(psf0,2)/2 )+1;
    
    %% Take the center partition for an easier display.
    psf_line_0 = psf0( :,iCenterPartition) ; % extract center line 
    psf_line_1 = psf1( :,iCenterPartition) ; % extract center line 
    psf_line_2 = psf2( :,iCenterPartition) ; % extract center line 
    psf_line_3 = psf3( :,iCenterPartition) ; % extract center line 

    [~, index] = max(psf_line_0);
    interXVec = interXVec( :,iCenterPartition);
    shortXVec = interXVec - interXVec(index);
    
    ylimitUP = max([psf_line_0;psf_line_1;psf_line_2], [],'all') *1.1;
    ylimitDOWN = min([psf_line_0;psf_line_1;psf_line_2], [],'all') *1.2;
    [xPos1, xPos2] = CR_getFWHM(shortXVec, psf_line_2);

    figure;
    plot( shortXVec, psf_line_0, 'LineWidth', 3)
    hold on
    plot( shortXVec, psf_line_1, 'LineWidth', 3)
    plot( shortXVec, psf_line_2, 'LineWidth', 3)
    xlabel('Voxel Index')
    ylabel('PSF Amplitude (a.u.)')
    title(strcat('Turbofactor = ',num2str(turboFactor(z))))
    ax = gca;    ax.FontSize = 16; 
    xline(xPos1,'-.') %
    xline(xPos2,'-.')
    xlim([-5 5])
    ylim([ylimitDOWN ylimitUP])
    xticks(-5:1:5)
    set(gcf,'position',[100,100,600,500])
    legend( 'Single', 'Dual', 'ihMTR')

    if SaveImgs; saveas(gcf,strcat(SavDir,'PSF_slice_',num2str(Params.TurboFactor),'turbofactor.png')); end

    
    % save the ihMT one to look at between turbo factors after!
    psf_single_values(z,:) = psf_line_0;
    psf_dual_values(z,:) = psf_line_1;
    psf_ihMT_values(z,:) = psf_line_2;
    psf_Ref_values(z,:) = psf_line_3;

pause(3);
close;

end % end of turbofactor loop



%% With all of them run, generate a plot with the different ihMTs and normalize to the greatest values

maxIHMT = max( psf_ihMT_values,[],'all');

% normalize each ihMT to this
temp = psf_ihMT_values;
temp = temp./maxIHMT;

% lowest turbofactor should have the narrowest FWHM, use that as reference.
[xPos1, xPos2] = CR_getFWHM(shortXVec', psf_ihMT_values(end,:));
% worst case scenario
[xPos3, xPos4] = CR_getFWHM(shortXVec', psf_ihMT_values(1,:));

figure;
clrVal = [0 0 0.9];
 plot(shortXVec, temp(1,:), 'LineWidth', 3, 'Color',clrVal)
hold on
for i = 2:size(temp,1)
    if i == length(turboFactor)
        clrVal = [1 0.9 0.1];
    else
        %clrVal = clrVal + [0.10 0.15 0.01];
        clrVal(1) = clrVal(1) +0.1;
        clrVal(2) = (1-clrVal(2))/2 + clrVal(2);
        clrVal(3) = clrVal(3) + 0.015;
    end
    plot(shortXVec, temp(i,:), 'LineWidth', 3, 'Color',clrVal)
end
xlabel('Voxel Index')
ylabel('PSF Amplitude (normalized)')
%title('Normalized ihMT PSF ')
ax = gca;    ax.FontSize = 20; 
    xline([xPos1 xPos2],'--','Color',[0 0.7 0],'LineWidth',2) %
    xline([xPos3 xPos4],'--','Color',[0.7 0 0],'LineWidth',2) %
  xlim([-5 5])
  xticks(-5:1:5)
set(gcf,'position',[100,100,800,600])
title('ihMTR PSF')
legend( 'TF = 200','TF = 160','TF = 120','TF = 80','TF = 48','TF = 10') 
if SaveImgs; saveas(gcf,strcat(SavDir,'ihMTR_normalized_ACROSS_turbofactors.png')); end


%% Repeat for Single

maxSingle = max( psf_single_values,[],'all');

% normalize each ihMT to this
temp = psf_single_values;
temp = temp./maxSingle;

% lowest turbofactor should have the narrowest FWHM, use that as reference.
[xPos1, xPos2] = CR_getFWHM(shortXVec', psf_single_values(end,:));
% worst case scenario
[xPos3, xPos4] = CR_getFWHM(shortXVec', psf_single_values(1,:));

figure;
clrVal = [0 0 0.9];
 plot(shortXVec, temp(1,:), 'LineWidth', 3, 'Color',clrVal)
hold on
for i = 2:size(temp,1)
    if i == length(turboFactor)
        clrVal = [1 0.9 0.1];
    else
        %clrVal = clrVal + [0.10 0.15 0.01];
        clrVal(1) = clrVal(1) +0.1;
        clrVal(2) = (1-clrVal(2))/2 + clrVal(2);
        clrVal(3) = clrVal(3) + 0.015;
    end
    plot(shortXVec, temp(i,:), 'LineWidth', 3, 'Color',clrVal)
end
xlabel('Voxel Index')
ylabel('PSF Amplitude (normalized)')
%title('Normalized ihMT PSF ')
ax = gca;    ax.FontSize = 20; 
    xline([xPos1 xPos2],'--','Color',[0 0.7 0],'LineWidth',2) %
    xline([xPos3 xPos4],'--','Color',[0.7 0 0],'LineWidth',2) %
  xlim([-5 5])
  xticks(-5:1:5)
set(gcf,'position',[100,100,800,600])
legend( 'TF = 200','TF = 160','TF = 120','TF = 80','TF = 48','TF = 10') 
title('Single PSF')
if SaveImgs; saveas(gcf,strcat(SavDir,'Single_normalized_ACROSS_turbofactors.png')); end

%% Repeat for Dual

maxDual = max( psf_dual_values,[],'all');

% normalize each ihMT to this
temp = psf_dual_values;
temp = temp./maxDual;

% lowest turbofactor should have the narrowest FWHM, use that as reference.
[xPos1, xPos2] = CR_getFWHM(shortXVec', psf_dual_values(end,:));
% worst case scenario
[xPos3, xPos4] = CR_getFWHM(shortXVec', psf_dual_values(1,:));

figure;
clrVal = [0 0 0.9];
 plot(shortXVec, temp(1,:), 'LineWidth', 3, 'Color',clrVal)
hold on
for i = 2:size(temp,1)
    if i == length(turboFactor)
        clrVal = [1 0.9 0.1];
    else
        %clrVal = clrVal + [0.10 0.15 0.01];
        clrVal(1) = clrVal(1) +0.1;
        clrVal(2) = (1-clrVal(2))/2 + clrVal(2);
        clrVal(3) = clrVal(3) + 0.015;
    end
    plot(shortXVec, temp(i,:), 'LineWidth', 3, 'Color',clrVal)
end
xlabel('Voxel Index')
ylabel('PSF Amplitude (normalized)')
%title('Normalized ihMT PSF ')
ax = gca;    ax.FontSize = 20; 
    xline([xPos1 xPos2],'--','Color',[0 0.7 0],'LineWidth',2) %
    xline([xPos3 xPos4],'--','Color',[0.7 0 0],'LineWidth',2) %
  xlim([-5 5])
  xticks(-5:1:5)
set(gcf,'position',[100,100,800,600])
legend( 'TF = 200','TF = 160','TF = 120','TF = 80','TF = 48','TF = 10') 
title('Dual PSF')
if SaveImgs; saveas(gcf,strcat(SavDir,'Dual_normalized_ACROSS_turbofactors.png')); end





%% Want to better visualize the comparison between SNR and PSF resolution. Plot SNR (y) vs FWHM (x)
v = size(psf_ihMT_values,1);
% First pull out FWHM.
x1 = zeros(2,v);
for i = 1:v
    [x1(1,i), x1(2,i)] = CR_getFWHM(shortXVec', psf_ihMT_values(i,:));
end
FWHM = x1(2,:) - x1(1,:); 

% The SNR values will come from my simulation results. Just list them here:
SNR = [ 7.209, 5.8235, 4.225,4.002, 3.595, 2.2965  ];
% Lets add a linear line that corresponds to the linear increase in SNR
% from voxel size alone. 
m = SNR(end)/FWHM(end);
b = 0;

% Make the plot
figure;
hold on
clrVal = [0 0 0.9];
ylim([0 max(SNR)*1.11])
xlim([0.5 3])
for i = 1:v
    if i == length(turboFactor)
        clrVal = [1 0.9 0.1];
        scatter(FWHM(i), SNR(i), 125,clrVal,'filled')
    else
        scatter(FWHM(i), SNR(i), 125,clrVal,'filled')
        %clrVal = clrVal + [0.10 0.15 0.01];
        clrVal(1) = clrVal(1) +0.1;
        clrVal(2) = (1-clrVal(2))/2 + clrVal(2);
        clrVal(3) = clrVal(3) + 0.015;
    end
    %scatter(FWHM(i), SNR(i), 125,clrVal,'filled')
end
refline(m,b)
scatter(FWHM(i), SNR(i), 125,clrVal,'filled')
ax = gca;    ax.FontSize = 20; 
set(gcf,'position',[100,100,800,600])
xlabel('FWHM (mm)')
ylabel('ihMTR SNR')
xline(1,'--',{'Nominal','Resolution'},'Color',[0 0.0 0],'LineWidth',2,'FontSize',16) %
legend( 'TF = 200','TF = 160','TF = 120','TF = 80','TF = 48','TF = 10','Location','southeast') 
%hold off

if SaveImgs; saveas(gcf,strcat(SavDir,'SNR_vs_FWHM.png')); end









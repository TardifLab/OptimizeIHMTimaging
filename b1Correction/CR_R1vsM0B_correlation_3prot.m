%% Generate the M0B mapping to R1 from simulation results and acquired data

addpath(genpath('Directory/niak-master' ))
addpath(genpath('Directory/NeuroImagingMatlab'))
addpath(genpath( 'Directory/MP2RAGE'   ))

OutputDir = 'Directory\b1Correction\outputs';
%% Load images:


DATADIR = 'Image/Directory';

%image names:
% in the order of dual, hfa, neg, lfa, pos
mtw_fn = {'dual_reg.mnc' 'pos_reg.mnc' 'neg_reg.mnc' 'noMT_reg.mnc' ...                           % ihMT
          'dual_boost_reg.mnc' 'pos_boost_reg.mnc' 'neg_boost_reg.mnc' 'noMT_boost_reg.mnc' ...   % boost v1
          'dual_boost2_reg.mnc' 'pos_boost2_reg.mnc' 'neg_boost2_reg.mnc' 'noMT_boost2_reg.mnc' ... % boost v2
          'spMP2RAGE_inv1_reg.mnc.gz' 'spMP2RAGE_inv2_reg.mnc.gz' 'spMP2RAGE_UNI_reg.mnc.gz'...     % Sparse MP2RAGE
          };                                            

for i = 1:size(mtw_fn,2)
    fn = fullfile(DATADIR,mtw_fn{i});
    [hdr, img] = niak_read_vol(fn);
    comb_mtw(:,:,:,i) = img; %.img;
end


%% Load the mask
fn = fullfile(DATADIR,'itkmask.nii.gz');
[~, mask] = niak_read_vol(fn);
mask1 = permute(mask,[2 3 1]); % conversion between minc and nii reorients it

%% Some B1 issues so lets try and load that
[~, b1] = niak_read_vol(fullfile(DATADIR,'resampled_b1field.mnc')); 
b1 = double(b1);
b1 = limitHandler(b1, 0.5,1.4);
b1 = CR_imgaussfilt3_withMask(b1, mask1, 5); %light smoothing to the map

figure; imshow3Dfull(b1, [0.6 1.2],jet)

%% Run MP-PCA denoising
comb_mtw = double(comb_mtw);
all_PCAcorr = MPdenoising(comb_mtw);


%% separate the images then average the MTw ones
dual = all_PCAcorr(:,:,:,1);
pos  = all_PCAcorr(:,:,:,2);
neg = all_PCAcorr(:,:,:,3);
noMT = all_PCAcorr(:,:,:,4);

dual2 = all_PCAcorr(:,:,:,5);
pos2  = all_PCAcorr(:,:,:,6);
neg2 = all_PCAcorr(:,:,:,7);
noMT2 = all_PCAcorr(:,:,:,8);

dual3 = all_PCAcorr(:,:,:,9);
pos3  = all_PCAcorr(:,:,:,10);
neg3 = all_PCAcorr(:,:,:,11);
noMT3 = all_PCAcorr(:,:,:,12);

sp_mp2r_inv1  = all_PCAcorr(:,:,:,13);
sp_mp2r_inv2 = all_PCAcorr(:,:,:,14);
sp_mp2r_uni  = all_PCAcorr(:,:,:,15);


%% Now from the MP2RAGE:
  
MP2RAGE.B0 = 3;           % in Tesla
MP2RAGE.TR = 5;           % MP2RAGE TR in seconds
MP2RAGE.TRFLASH = 6.4e-3; % TR of the GRE readout
MP2RAGE.TIs = [0.94 2.83];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
MP2RAGE.NZslices = [ceil(175/2) floor(175/2)];%  should be two values, number of excitations before k-space center, and number after. [Slices Per Slab * [PartialFourierInSlice-0.5  0.5] ]
MP2RAGE.FlipDegrees = [4 5];% Flip angle of the two readouts in degrees


MP2RAGEimg.img = sp_mp2r_uni; % load_untouch_nii(MP2RAGE.filenameUNI);
MP2RAGEINV2img.img = sp_mp2r_inv2; % load_untouch_nii(MP2RAGE.filenameINV2);
B1.img = b1;
brain.img = mask1;

tic
[ spT1map, spMP2RAGEcorrected, spAppmap2] = CR_T1B1correctpackageTFL_withM0( B1, MP2RAGEimg, MP2RAGEINV2img, MP2RAGE, brain, 0.96);
toc

spT1_map = spT1map.img;
spApp_mp2 = spAppmap2.img ;

spT1_map = limitHandler(spT1_map);
spApp_mp2 = double(limitHandler(spApp_mp2));

figure; imshow3Dfull(spT1_map, [300 2500],jet)
figure; imshow3Dfull(spApp_mp2 , [00 6000])

hdr.file_name = fullfile(DATADIR,'matlab/sparseMP2RAGE_T1.mnc.gz'); niak_write_vol(hdr, spT1_map);
hdr.file_name = fullfile(DATADIR,'matlab/sparseMP2RAGE_M0.mnc.gz'); niak_write_vol(hdr, spApp_mp2);


%% Mask -> bet result touched up in itk, then threshold CSF and some dura

mask = mask1;
mask(spT1_map > 2500) = 0;
mask(spT1_map < 650) = 0;
mask(isnan(spT1_map)) = 0;
mask = bwareaopen(mask, 10000,6);
figure; imshow3Dfullseg(spT1_map, [300 2500],mask)


%% Compute MTsat ihMTsat

%% Protocol 1
echoSpacing = 7.66; % The echo spacing of the GRE readout
numExcitation = 10; 
TR = 120;
flipA = 5; % flip angle
DummyEcho = 2;

sat_dual1 = calcMTsatThruLookupTablewithDummyV3( dual, b1, spT1_map, mask1,spApp_mp2, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_pos1  = calcMTsatThruLookupTablewithDummyV3( pos, b1, spT1_map, mask1, spApp_mp2, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_neg1  = calcMTsatThruLookupTablewithDummyV3( neg, b1, spT1_map, mask1, spApp_mp2, echoSpacing, numExcitation, TR, flipA, DummyEcho);

% figure; imshow3Dfull(sat_dual1 , [0 0.06], jet); figure; imshow3Dfull(sat_pos1 , [0 0.05], jet);  figure; imshow3Dfull(sat_neg1 , [0 0.05], jet); 

%% Protocol 2
numExcitation = 202; 
TR = 3000;
flipA = 11; % flip angle

sat_dual2 = calcMTsatThruLookupTablewithDummyV3( dual2, b1, spT1_map, mask1,spApp_mp2, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_pos2  = calcMTsatThruLookupTablewithDummyV3( pos2, b1, spT1_map, mask1, spApp_mp2, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_neg2  = calcMTsatThruLookupTablewithDummyV3( neg2, b1, spT1_map, mask1, spApp_mp2, echoSpacing, numExcitation, TR, flipA, DummyEcho);

% figure; imshow3Dfull(sat_dual2 , [0 0.46], jet); figure; imshow3Dfull(sat_pos2 , [0 0.46], jet) ; figure; imshow3Dfull(sat_neg2 , [0 0.46], jet); 

%% Protocol 3
numExcitation = 82; 
TR = 1140;
flipA = 7; % flip angle

sat_dual3 = calcMTsatThruLookupTablewithDummyV3( dual3, b1, spT1_map, mask1,spApp_mp2, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_pos3  = calcMTsatThruLookupTablewithDummyV3( pos3, b1, spT1_map, mask1, spApp_mp2, echoSpacing, numExcitation, TR, flipA, DummyEcho);
sat_neg3  = calcMTsatThruLookupTablewithDummyV3( neg3, b1, spT1_map, mask1, spApp_mp2, echoSpacing, numExcitation, TR, flipA, DummyEcho);

% figure; imshow3Dfull(sat_dual3 , [0 0.36], jet);


%% With MTsat maps made, perform M0b mapping


% load in the fit results for VFA - Optimal
fitValues_S_8 = load(fullfile(OutputDir,'fitValues_S_8.mat'));
fitValues_S_8 = fitValues_S_8.fitValues;
fitValues_D_8 = load(fullfile(OutputDir,'fitValues_D_8.mat'));
fitValues_D_8 = fitValues_D_8.fitValues;


% % load in the fit results for VFA - Boosted approach
fitValues_S_200 = load(fullfile(OutputDir,'fitValues_S_200.mat'));
fitValues_S_200 = fitValues_S_200.fitValues;
fitValues_D_200 = load(fullfile(OutputDir,'fitValues_D_200.mat'));
fitValues_D_200 = fitValues_D_200.fitValues;

% 
% % load in the fit results for VFA - MTsat approach
fitValues_S_80 = load(fullfile(OutputDir,'fitValues_S_80.mat'));
fitValues_S_80 = fitValues_S_80.fitValues;
fitValues_D_80 = load(fullfile(OutputDir,'fitValues_D_80.mat'));
fitValues_D_80 = fitValues_D_80.fitValues;


% need to convert to 1/s from 1/ms -> ONLY USE MP2RAGE values, VFA are too
% far off.
R1_s = (1./spT1_map) *1000;

% initialize matrices
M0b_8_dual = zeros(size(sat_dual1));
M0b_8_pos = zeros(size(sat_dual1));
M0b_8_neg = zeros(size(sat_dual1));

M0b_200_dual = zeros(size(sat_dual1));
M0b_200_pos = zeros(size(sat_dual1));
M0b_200_neg = zeros(size(sat_dual1));

M0b_80_dual = zeros(size(sat_dual1));
M0b_80_pos = zeros(size(sat_dual1));
M0b_80_neg = zeros(size(sat_dual1));


%% SPEED IT UP BY DOING A FEW AXIAL SLICES
axialStart = 126; % 65
axialStop = axialStart+3;%115;
% check
 figure; imshow3Dfull(sat_dual1(:,axialStart:axialStop,:) , [0 0.06], jet)

b1_1 = 11.4;
b1_2 = 13.3;
b1_3 = 11.6;

tic %  
for i = 1:size(sat_dual1,1) % went to 149
    
    for j = axialStart:axialStop % 1:size(sat_dual,2) % for axial slices
        for k =  1:size(sat_dual1,3) % sagital slices  65
            
            if mask(i,j,k) > 0 %&& dual_s(i,j,k,3) > 0
                
                 
                 [M0b_8_dual(i,j,k), ~,  ~]  = CR_fit_M0b_v2( b1_1*b1(i,j,k), R1_s(i,j,k), sat_dual1(i,j,k), fitValues_D_8);
                 [M0b_8_pos(i,j,k),  ~,  ~]  = CR_fit_M0b_v2( b1_1*b1(i,j,k), R1_s(i,j,k), sat_pos1(i,j,k), fitValues_S_8);               
                 [M0b_8_neg(i,j,k),  ~,  ~]  = CR_fit_M0b_v2( b1_1*b1(i,j,k), R1_s(i,j,k), sat_neg1(i,j,k), fitValues_S_8);
                 
                 % TF = 200 protocol
                 [M0b_200_dual(i,j,k), ~, ~] = CR_fit_M0b_v2( b1_2*b1(i,j,k), R1_s(i,j,k), sat_dual2(i,j,k),fitValues_D_200);
                 [M0b_200_pos(i,j,k),  ~, ~] = CR_fit_M0b_v2( b1_2*b1(i,j,k), R1_s(i,j,k), sat_pos2(i,j,k),fitValues_S_200);
                 [M0b_200_neg(i,j,k),  ~, ~] = CR_fit_M0b_v2( b1_2*b1(i,j,k), R1_s(i,j,k), sat_neg2(i,j,k),fitValues_S_200);
                 
                 % TF = 80  protocol
                 [M0b_80_dual(i,j,k), ~, ~] = CR_fit_M0b_v2( b1_3*b1(i,j,k), R1_s(i,j,k), sat_dual3(i,j,k),fitValues_D_80);
                 [M0b_80_pos(i,j,k),  ~, ~] = CR_fit_M0b_v2( b1_3*b1(i,j,k), R1_s(i,j,k), sat_pos3(i,j,k),fitValues_S_80);
                 [M0b_80_neg(i,j,k),  ~, ~] = CR_fit_M0b_v2( b1_3*b1(i,j,k), R1_s(i,j,k), sat_neg3(i,j,k),fitValues_S_80);     
                 
            end
        end
    end
    disp(i)
end
toc %% this took 30hours for 1mm isotropic full brain dataset. * was running fitting in another matlab
    % instance, so could be easily sped up running on its own and/or adding
    % the parfor loop. 


 figure; imshow3Dfull(M0b_8_pos, [0 0.15],jet)
  figure; imshow3Dfull(M0b_80_pos, [0 0.15],jet)
 figure; imshow3Dfull(M0b_200_pos, [0 0.15],jet)
    

figure; imshow3Dfull(M0b_200_dual, [0 0.15],jet)
figure; imshow3Dfull(M0b_8_dual, [0 0.15],jet)
figure; imshow3Dfull(M0b_80_dual, [0 0.15],jet)


% export
mkdir(fullfile(OutputDir,'processing'))
hdr.file_name = fullfile(OutputDir,'processing/M0b_8_dual.mnc.gz'); niak_write_vol(hdr,M0b_8_dual);
hdr.file_name = fullfile(OutputDir,'processing/M0b_8_pos.mnc.gz'); niak_write_vol(hdr,M0b_8_pos);
hdr.file_name = fullfile(OutputDir,'processing/M0b_8_neg.mnc.gz'); niak_write_vol(hdr,M0b_8_neg);

hdr.file_name = fullfile(OutputDir,'processing/M0b_200_dual.mnc.gz'); niak_write_vol(hdr,M0b_200_dual);
hdr.file_name = fullfile(OutputDir,'processing/M0b_200_pos.mnc.gz'); niak_write_vol(hdr,M0b_200_pos);
hdr.file_name = fullfile(OutputDir,'processing/M0b_200_neg.mnc.gz'); niak_write_vol(hdr,M0b_200_neg);

hdr.file_name = fullfile(OutputDir,'processing/M0b_80_dual.mnc.gz'); niak_write_vol(hdr,M0b_80_dual);
hdr.file_name = fullfile(OutputDir,'processing/M0b_80_pos.mnc.gz'); niak_write_vol(hdr,M0b_80_pos);
hdr.file_name = fullfile(OutputDir,'processing/M0b_80_neg.mnc.gz'); niak_write_vol(hdr,M0b_80_neg);


%% With M0B maps made, correlate with R1 and update the fitValues file. 

% use this fake mask to get rid of dura. 
tempMask = mask;
tempMask = imerode(tempMask, strel('sphere',2));
figure; imshow3Dfullseg(M0b_200_dual, [0 0.15],tempMask)

mkdir(fullfile(OutputDir,'figures'));

% Optimized Approach
fitValues_D_8  = CR_generate_R1vsM0B_correlation( R1_s, M0b_8_dual, tempMask, fitValues_D_8, fullfile(OutputDir,'figures/R1vsM0b_8_dual.png'), fullfile(OutputDir,'fitValues_D_8.mat'));
fitValues_SP_8 = CR_generate_R1vsM0B_correlation( R1_s, M0b_8_pos, tempMask, fitValues_S_8, fullfile(OutputDir,'figures/R1vsM0b_8_pos.png'), fullfile(OutputDir,'fitValues_SP_8.mat'));
fitValues_SN_8 = CR_generate_R1vsM0B_correlation( R1_s, M0b_8_neg, tempMask, fitValues_S_8, fullfile(OutputDir,'figures/R1vsM0b_8_neg.png'), fullfile(OutputDir,'fitValues_SN_8.mat'));


% Boosted Approach
fitValues_D_200  = CR_generate_R1vsM0B_correlation( R1_s, M0b_200_dual, tempMask, fitValues_D_200, fullfile(OutputDir,'figures/R1vsM0b_200_dual.png'), fullfile(OutputDir,'fitValues_D_200.mat'));
fitValues_SP_200 = CR_generate_R1vsM0B_correlation( R1_s, M0b_200_pos, tempMask, fitValues_S_200, fullfile(OutputDir,'figures/R1vsM0b_200_pos.png'), fullfile(OutputDir,'fitValues_SP_200.mat'));
fitValues_SN_200 = CR_generate_R1vsM0B_correlation( R1_s, M0b_200_neg, tempMask, fitValues_S_200, fullfile(OutputDir,'figures/R1vsM0b_200_neg.png'), fullfile(OutputDir,'fitValues_SN_200.mat'));

% MTsat Approach
fitValues_D_80  = CR_generate_R1vsM0B_correlation( R1_s, M0b_80_dual, tempMask, fitValues_D_80, fullfile(OutputDir,'figures/R1vsM0b_80_dual.png'), fullfile(OutputDir,'fitValues_D_80.mat'));
fitValues_SP_80 = CR_generate_R1vsM0B_correlation( R1_s, M0b_80_pos, tempMask, fitValues_S_80, fullfile(OutputDir,'figures/R1vsM0b_80_pos.png'), fullfile(OutputDir,'fitValues_SP_80.mat'));
fitValues_SN_80 = CR_generate_R1vsM0B_correlation( R1_s, M0b_80_neg, tempMask, fitValues_S_80, fullfile(OutputDir,'figures/R1vsM0b_80_neg.png'), fullfile(OutputDir,'fitValues_SN_80.mat'));



%% Now use these results to B1 correct the data:
OutputDir = DATADIR;


b1_1 = 11.4;
b1_2 = 13.3;
b1_3 = 11.6;

corr_prot1_d = MTsat_B1corr_factor_map(b1, R1_s, b1_1, fitValues_D_8);
corr_prot1_p = MTsat_B1corr_factor_map(b1, R1_s, b1_1, fitValues_SP_8);
corr_prot1_n = MTsat_B1corr_factor_map(b1, R1_s, b1_1, fitValues_SN_8);

corr_prot2_d = MTsat_B1corr_factor_map(b1, R1_s, b1_2, fitValues_D_200);
corr_prot2_p = MTsat_B1corr_factor_map(b1, R1_s, b1_2, fitValues_SP_200);
corr_prot2_n = MTsat_B1corr_factor_map(b1, R1_s, b1_2, fitValues_SN_200);

corr_prot3_d = MTsat_B1corr_factor_map(b1, R1_s, b1_3, fitValues_D_80);
corr_prot3_p = MTsat_B1corr_factor_map(b1, R1_s, b1_3, fitValues_SP_80);
corr_prot3_n = MTsat_B1corr_factor_map(b1, R1_s, b1_3, fitValues_SN_80);


% Part 2, apply correction map
sat_dual1_c = (sat_dual1 + sat_dual1.* corr_prot1_d) .* mask1;
sat_pos1_c  = (sat_pos1 + sat_pos1.* corr_prot1_p) .* mask1;
sat_neg1_c  = (sat_neg1 + sat_neg1.* corr_prot1_n) .* mask1;
ihmt1_c      = sat_dual1_c - (sat_pos1_c + sat_neg1_c)/2;

sat_dual2_c = (sat_dual2 + sat_dual2.* corr_prot2_d) .* mask1;
sat_pos2_c  = (sat_pos2 + sat_pos2.* corr_prot2_p) .* mask1;
sat_neg2_c  = (sat_neg2 + sat_neg2.* corr_prot2_n) .* mask1;
ihmt2_c      = sat_dual2_c - (sat_pos2_c + sat_neg2_c)/2;

sat_dual3_c = (sat_dual3 + sat_dual3.* corr_prot3_d) .* mask1;
sat_pos3_c  = (sat_pos3 + sat_pos3.* corr_prot3_p) .* mask1;
sat_neg3_c  = (sat_neg3 + sat_neg3.* corr_prot3_n) .* mask1;
ihmt3_c      = sat_dual3_c - (sat_pos3_c + sat_neg3_c)/2;


ihmt1_c = double(limitHandler(ihmt1_c,0, 0.05));
ihmt2_c = double(limitHandler(ihmt2_c,0, 0.15));
ihmt3_c = double(limitHandler(ihmt3_c,0, 0.10));


ihmt1_c( ihmt1_c >= 0.05) = 0;
ihmt2_c( ihmt2_c >= 0.15) = 0;
ihmt3_c( ihmt3_c >= 0.1) = 0;

%% View results
figure; imshow3Dfull(sat_dual1_c , [0 0.06], jet); 
figure; imshow3Dfull(sat_pos1_c , [0 0.06], jet)
figure; imshow3Dfull(sat_neg1_c , [0 0.06], jet); 

figure; imshow3Dfull(sat_dual2_c , [0 0.5], jet); 
figure; imshow3Dfull(sat_pos2_c , [0 0.5], jet)
figure; imshow3Dfull(sat_neg2_c , [0 0.5], jet); 

figure; imshow3Dfull(sat_dual3_c , [0 0.35], jet); 
figure; imshow3Dfull(sat_pos3_c , [0 0.35], jet)
figure; imshow3Dfull(sat_neg3_c , [0 0.35], jet); 

figure; imshow3Dfull(ihmt1_c , [0 0.06], jet)
figure; imshow3Dfull(ihmt2_c , [0 0.15], jet)
figure; imshow3Dfull(ihmt3_c , [0 0.08], jet)

%% Other things, save if you want

hdr.file_name = strcat(DATADIR,'matlab/MTsat_dual_1.mnc.gz'); niak_write_vol(hdr, sat_dual1_c);
hdr.file_name = strcat(DATADIR,'matlab/MTsat_pos_1.mnc.gz'); niak_write_vol(hdr, sat_pos1_c);
hdr.file_name = strcat(DATADIR,'matlab/MTsat_neg_1.mnc.gz'); niak_write_vol(hdr, sat_neg1_c);
hdr.file_name = strcat(DATADIR,'matlab/ihMTsat_1.mnc.gz'); niak_write_vol(hdr, ihmt1_c);

hdr.file_name = strcat(DATADIR,'matlab/MTsat_dual_2.mnc.gz'); niak_write_vol(hdr, sat_dual2_c);
hdr.file_name = strcat(DATADIR,'matlab/MTsat_pos_2.mnc.gz'); niak_write_vol(hdr, sat_pos2_c);
hdr.file_name = strcat(DATADIR,'matlab/MTsat_neg_2.mnc.gz'); niak_write_vol(hdr, sat_neg2_c);
hdr.file_name = strcat(DATADIR,'matlab/ihMTsat_2.mnc.gz'); niak_write_vol(hdr, ihmt2_c);

hdr.file_name = strcat(DATADIR,'matlab/MTsat_dual_3.mnc.gz'); niak_write_vol(hdr, sat_dual3_c);
hdr.file_name = strcat(DATADIR,'matlab/MTsat_pos_3.mnc.gz'); niak_write_vol(hdr, sat_pos3_c);
hdr.file_name = strcat(DATADIR,'matlab/MTsat_neg_3.mnc.gz'); niak_write_vol(hdr, sat_neg3_c);
hdr.file_name = strcat(DATADIR,'matlab/ihMTsat_3.mnc.gz'); niak_write_vol(hdr, ihmt3_c);


hdr.file_name = strcat(DATADIR,'matlab/b1.mnc.gz'); niak_write_vol(hdr, b1);



%% 



ihmtSlice1 = squeeze( ihmt1_c(:,126,:));
ihmtSlice2 = squeeze( ihmt2_c(:,126,:));
ihmtSlice3 = squeeze( ihmt3_c(:,126,:));

figure; imagesc(ihmtSlice1); axis image;
colormap(gray)
caxis([0 0.03])
hold on
line([0,175], [165,166], 'Color', 'r');


ihmtProf1 = ihmtSlice1(165,:);
ihmtProf2 = ihmtSlice2(165,:);
ihmtProf3 = ihmtSlice3(165,:);

% normalize:
ihmtProf1 = ihmtProf1/ max(ihmtProf1);
ihmtProf2 = ihmtProf2/ max(ihmtProf2);
ihmtProf3 = ihmtProf3/ max(ihmtProf3);


figure
plot(ihmtProf1,'LineWidth',2); 
hold on;
plot(ihmtProf3,'LineWidth',2); 
plot(ihmtProf2,'LineWidth',2); 
title('Line Profile (L-R)');
legend(' TF = 8', ' TF = 80', ' TF = 200', 'Location','northwest');
xlim([15, 160])
hold off
xlabel('Voxel Index');
ylabel('Relative ihMT_{sat}')
ax = gca; ax.FontSize = 20; 
set(gcf,'position',[10,400,1000,600])   





















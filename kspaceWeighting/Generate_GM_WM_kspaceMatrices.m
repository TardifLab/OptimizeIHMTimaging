%% Generate_GM_WM_kspaceMatrices
% This is mainly needed for centric-out encoding on Siemens.

DATADIR = 'kspaceWeighting/Atlas_reference/mni_icbm152_nlin_sym_09a_minc2/';

% Load the maps
mtw_fn = {'mni_icbm152_csf_tal_nlin_sym_09a.mnc';'mni_icbm152_gm_tal_nlin_sym_09a.mnc'; 'mni_icbm152_wm_tal_nlin_sym_09a.mnc';'DeepStructureMask.mnc'};
                                                
% load images
%comb_mtw = zeros(224, 256, 176, 18);

for i = 1:length(mtw_fn)
    fn = strcat(DATADIR,mtw_fn{i});
    %[hdr, img] = niak_read_vol(fn);
    [hdr, img] = minc_read(fn);
    comb_mtw(:,:,:,i) = img; %.img;
end

figure; imshow3Dfull(comb_mtw(:,:,:,1)  )

%% Pre-emptively change the matrix size to match what I collected (176x216x192)

comb_mtw2 = zeros( 176,216,192,length(mtw_fn) );
comb_mtw2(:,:,2:190,:) = comb_mtw( 11:186 , 9:224,:, :);

[x,y,z] = size(comb_mtw2(:,:,:,1));

comb_mtw = comb_mtw2;
clearvars comb_mtw2;

switch Params.Orientation

    case  'Axial'
        comb_mtw2 = comb_mtw; % already in this 
    case 'Sagittal'
        comb_mtw2 = permute(comb_mtw, [2,3,1,4]);
    case 'Coronal'
        comb_mtw2 = permute(comb_mtw, [3,1,2,4]);
    otherwise 
        disp( 'Set Params.Orientation to either Axial, Sagittal, or Coronal')
end

[x,y,z] = size(comb_mtw2(:,:,:,1));

comb_mtw = comb_mtw2;
clearvars comb_mtw2;


%% Brain tissue:
brain_m = zeros(x,y,z);
brain_m(comb_mtw(:,:,:,2) >=0.05) = 1; % add extra, then remove
brain_m(comb_mtw(:,:,:,3) >=0.05) = 1;
fft_brain_m = fftshift(fftn(brain_m));
figure; imshow3Dfull( brain_m )

%% Save results.

save( strcat(DATADIR,'GM_seg_MNI_152_image.mat'), 'brain_m')
save( strcat(DATADIR,'GM_seg_MNI_152_kspace.mat'), 'fft_brain_m')







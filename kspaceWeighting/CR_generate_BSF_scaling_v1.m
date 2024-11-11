function sig = CR_generate_BSF_scaling_v1( inputMag, Params, outputSamplingTable , brain_m, fft_brain_m)

% Over long excitation trains, you have a change in the magnetization.
% if you have a map that is being calculated off of the first echo in the
% train, it is assuming the signal == across the train.

% gm_m is the grey matter segmentation that is the same reconstructred size as the sampling table
% fft_gm_m is the k-space values of the above (complex!). It is expected
% that fftshift has already been run on this value (middle is the
% brightest)


%% Fill table with echo train values
outKspace_s = CR_fillKspaceSamplingTable_v2( inputMag, outputSamplingTable, Params);

%% Interpolate missing grappa lines
outKspace_s = CR_interpolateMissingGrappaLines( outKspace_s);

%% Assume constant values in readout direction 
sim3d_m = repmat(outKspace_s, [1,1,Params.Slices]);

%% Scale the k-space by brain weighting to get 'Brain-spread-function'
bsf = sim3d_m .* fft_brain_m;

b_v = abs(ifftn(ifftshift(bsf)));

sig = mean( b_v( brain_m>0 ));



% Code for viewing/debugging
% figure; imagesc( outKspace_s);
% axis image; colorbar;

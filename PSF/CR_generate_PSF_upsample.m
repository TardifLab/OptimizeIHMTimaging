function [ Interpolated_PSF, interXVec, interYVec] = CR_generate_PSF_upsample( inputKspace )

% Takes an input k-space data and generates the corresponding PSF.
% The PSF is upsampled for your viewing pleasure

% Get size limits
[x, y] = size(inputKspace);
xVec =  -1*(x/2)+1: (x/2);     % X vectors to center on 0
yVec =  -1*(y/2)+1: (y/2);     % X vectors to center on 0

% Calculate psf
psf0 = abs(ifftshift(ifft2(inputKspace))); % working with magnitude data, so use this
%figure; imagesc(psf0);colorbar; axis image;

% Interpolate the results for viewing...
IncrementVal = 0.01;
interXVec = -1*(x/2)+1 :IncrementVal: (x/2)+1; 
interYVec = -1*(y/2)+1 :IncrementVal: (y/2)+1; 

[xVecM, yVecM] = ndgrid(xVec,yVec);
[interXVec,interYVec] = ndgrid(interXVec,interYVec);
 
Interpolated_PSF = interpn(xVecM, yVecM, psf0, interXVec, interYVec, 'spline');



% figure; plot(Interpolated_PSF)
function [outputSamplingTable, measuredElem, iSegments] = Step1_calculateKspaceSampling_v3 (Params)

%% Make an ordering scheme for elliptical-centric out. %

% Eventually turn this into a function that takes a Params input with the
% below values:

% iNumLines        = 64; % divide by 2 for grappa =2 factor
% iNumPartitions   = 16;
% iCenterLine      = floor(iNumLines/2 );
% iCenterPartition = floor(iNumPartitions/2);
% iSegments        = 8;
% iNumElements = double(iNumLines*iNumPartitions);

if Params.Grappa
    iNumLines        = Params.NumLines/Params.AccelerationFactor +Params.ReferenceLines; 
else
    iNumLines        = Params.NumLines;
end

% Added for consistency with previous code
if ~isfield(Params, 'ellipMask')
    Params.ellipMask = 1;
end

iNumPartitions   = Params.NumPartitions; 
iCenterLine      = floor(iNumLines/2 );
iCenterPartition = floor(iNumPartitions/2);
iTurboFactor        = Params.TurboFactor; 
iNumElements = double(iNumLines*iNumPartitions);

%% first generate the segmentation and order generating functions
[AngleOrderFunc, RadiusOrderFunc] = GeneratePieSegFunOrdFun(iNumLines, iNumPartitions, iCenterLine, iCenterPartition, iNumElements);

%% second, use the above functoins to generate the order table:

if Params.ellipMask
    A = ellipMask([iCenterLine*1.02 iCenterPartition*1.02], [iNumLines iNumPartitions], [iCenterLine iCenterPartition]) ;
    mask = A(:);
else
    mask = ones(size(RadiusOrderFunc));   
end

[outputSamplingTable, measuredElem, iSegments] = calcGeneratingFunctionTables_v3(AngleOrderFunc, RadiusOrderFunc,iNumElements,iTurboFactor,mask);

outputSamplingTable(mask == 0) = 0;


%% Can check results using

% viewResult = reshape(outputSamplingTable,Params.NumLines/Params.AccelerationFactor +Params.ReferenceLines,Params.NumPartitions);
% figure;
% imagesc(viewResult)
% axis image
% %caxis([0 15])
% colormap(jet)


% viewResult = reshape( raiOrderGeneratingFunction,Params.NumLines/Params.AccelerationFactor +Params.ReferenceLines,iNumPartitions);
% figure;
% imagesc(viewResult)
% axis image
% 
% viewResult = reshape(outputSamplingTable,iNumLines,iNumPartitions);
% figure;
% imagesc(viewResult)
% axis image

% viewResult = reshape( mask,Params.NumLines/Params.AccelerationFactor +Params.ReferenceLines,iNumPartitions);
% figure;
% imagesc(viewResult)
% axis image
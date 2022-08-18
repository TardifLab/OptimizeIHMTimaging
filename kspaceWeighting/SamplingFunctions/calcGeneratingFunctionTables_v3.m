function [FinalMat, m_iNumberOfMeasuredElements, iSegments] = calcGeneratingFunctionTables_v3(AngleOrderFunc, RadiusOrderFunc, iNumElements, iTurboFactor, mask) 



maskedElements = sum(~mask); 
    
m_iNumberOfMeasuredElements = iNumElements - maskedElements;

% Sort out the segmentation based on remaining lines
iSegments = ceil(m_iNumberOfMeasuredElements ./  iTurboFactor);


%% Divide up kspace by angle
% map can remain unmasked for now. 
angleSeg = -1*ones(size(AngleOrderFunc));
AngleOrderMask = AngleOrderFunc;
AngleOrderMask(mask == 0) = 10000;

[~, sortIdx] = sort(AngleOrderMask);

id1 = 1;
id2 = iTurboFactor;

for i = 1:iSegments
    angleSeg(sortIdx(id1:id2)) = i;
    id1 = id2;
    id2 = id2+ iTurboFactor;
end

angleSeg(mask == 0) = 0; % remove extras at end of last excitation train

%% Fill by Radius

FinalMat = zeros(size(RadiusOrderFunc));
idx = 0;
% for each segment
for i = 1:iSegments

    % Easiest way is to use find on whole matrix. So set values we don't
    % want included very high.
    temp = RadiusOrderFunc;
    temp(angleSeg ~= i) = 1e5;
    [~, sortIdx] = sort(temp);

    % fill them by their order in the sort. 
    for j = 1:iTurboFactor
        if (idx == m_iNumberOfMeasuredElements)
            break;
        end

        FinalMat(sortIdx(j)) = j;
        idx = idx+1;
    end
end

% 
% viewResult = reshape(angleSeg,iNumLines,iNumPartitions);
% figure; imagesc(viewResult); axis image; colorbar; colormap('jet'); %ax = gca; ax.CLim = [0 9];
% 
% viewResult = reshape(RadiusOrderFunc,iNumLines,iNumPartitions);
% figure; imagesc(viewResult); axis image; colorbar; colormap('jet'); %ax = gca; ax.CLim = [0 9];
% 
% viewResult = reshape(FinalMat,iNumLines,iNumPartitions);
% figure; imagesc(viewResult); axis image; colorbar; colormap('jet'); %ax = gca; ax.CLim = [0 9];







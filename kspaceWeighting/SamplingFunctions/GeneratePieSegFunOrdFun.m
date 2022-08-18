function [raiSegmentationGeneratingFunction, raiOrderGeneratingFunction] = GeneratePieSegFunOrdFun(iNumLines,iNumPartitions,iCenterLine,iCenterPartition, iNumElements) 

%% Generate raiSegmentationGeneratingFunction, raiOrderGeneratingFunction
    
    raiSegmentationGeneratingFunction = zeros(iNumElements,1);
    raiOrderGeneratingFunction = zeros(iNumElements,1);
    iCounter = 1;
    
    for iPartitionId = 1:iNumPartitions 
    
        dPosX = double(iPartitionId - iCenterPartition);
        dSquaredPosX = double(dPosX * dPosX);
    
        for iLineId = 1:iNumLines
            dPosY        = double(iLineId - iCenterLine);
            dSquaredPosY = dPosY * dPosY;
    
            % calculate radius
            dRadius = sqrt (dSquaredPosX + dSquaredPosY);
    
            % calculate angle
            dAngle = atan2 (dPosY, dPosX);
    
            % get angle in [0, 2*pi]
            if (dAngle < 0)
                dAngle =dAngle + 2 * pi;
            end
    
            raiSegmentationGeneratingFunction (iCounter) = (dAngle);
            raiOrderGeneratingFunction        (iCounter) = (dRadius);
            iCounter = iCounter + 1;
        end
    end
    
    
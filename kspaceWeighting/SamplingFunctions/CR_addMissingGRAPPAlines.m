function outputTable = CR_addMissingGRAPPAlines( inputTable, Params)

% Params.NumLines = 256;
% Params.NumPartitions = 176; 
% Params.Grappa = 1;
% Params.ReferenceLines = 32;
% Params.AccelerationFactor = 2;
% inputTable = repmat( [1:(Params.NumLines/Params.AccelerationFactor +Params.ReferenceLines)]', 1, Params.NumPartitions);

inputTable = reshape( inputTable,Params.NumLines/Params.AccelerationFactor +Params.ReferenceLines ,Params.NumPartitions);
outputTable = zeros( Params.NumLines, Params.NumPartitions);

% Num Lines from center outwards without skipping == Reference Lines
% because you are already acquiring the other lines, so divide 2 cancels

if Params.Grappa
    
    midLineUpOut = Params.NumLines /2 +1;
    midLineUpIn = (Params.NumLines/Params.AccelerationFactor +Params.ReferenceLines) /2 +1;

    midLineDownOut = Params.NumLines /2;
    midLineDownIn = (Params.NumLines/Params.AccelerationFactor +Params.ReferenceLines) /2;
 
    % Move down from middle
    idx = 0;
    for i = 0:(Params.NumLines/Params.AccelerationFactor +Params.ReferenceLines) /2  -1

        outputTable(midLineDownOut - idx,:) = inputTable(midLineDownIn - i,:); 
        outputTable(midLineUpOut + idx,:) = inputTable(midLineUpIn + i,:); 

        % reference line section, increase index by 1
        % Note, don't divide by acceleration factor here because reference lines fills this
        % area
        if i < Params.ReferenceLines
            idx = idx+1;
        else
            idx = idx+Params.AccelerationFactor;
        end
    end
else 
    outputTable = inputTable;
end

% figure;
% imagesc(outputTable)
% axis image
% colormap(jet)
% 
% figure;
% imagesc(inputTable)
% axis image
% colormap(jet)





































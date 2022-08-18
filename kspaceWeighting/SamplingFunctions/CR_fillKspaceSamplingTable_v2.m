function kSpaceFill = CR_fillKspaceSamplingTable_v2( inputMagTrain, outputSamplingTable, Params)
% v2 switched to only using radius, paired with function Step1_calculateKspaceSampling_v3

fillSpace = zeros(size(outputSamplingTable));

for j = 1:Params.TurboFactor

    fillSpace( outputSamplingTable == j ) = inputMagTrain(j);

end

if Params.Grappa
    kSpaceFill = reshape(fillSpace, Params.NumLines/Params.AccelerationFactor +Params.ReferenceLines, Params.NumPartitions );
    kSpaceFill = CR_addMissingGRAPPAlines( kSpaceFill, Params);
else
    kSpaceFill = reshape(fillSpace, Params.NumLines, Params.NumPartitions );
end



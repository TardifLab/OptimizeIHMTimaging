function Vq = CR_interpolateMissingGrappaLines( inputTable)

[x, y] = size(inputTable);

% Generate index for acquired lines
AcqIdx = max(inputTable,[],2) > 0;
xVec = find(AcqIdx~=0);

% remove zeroline
inputTable(AcqIdx==0,:) = [];

% Generate grid points for acquired data

[yg, xg] = ndgrid(xVec,1:y);

[yfinal_g, xfinal_g] = ndgrid(1:x,1:y);

Vq = interp2(xg, yg,inputTable,xfinal_g, yfinal_g,'nearest', 0); % extrapolation value is 0

Vq(Vq <0.0001) = 0;



% figure;
% imagesc(inputTable)
% axis image
% colormap(jet)
% 
% figure;
% imagesc(Vq)
% axis image
% colormap(jet)
% 
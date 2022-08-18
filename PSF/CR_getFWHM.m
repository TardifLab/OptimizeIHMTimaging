function [xPos1, xPos2] = CR_getFWHM(xVec, yVec)

yVec = yVec/(max(yVec)) *100; % normalize it.

% get max value
maxVal = max(yVec);
halfVal = maxVal /2;

% find index values where closest...
mid = yVec- halfVal;

% pick a threshold.
upT = 0;
lowT = 0;
xPos1 = [];
xPos2 = [];
idx = 0;

while (isempty(xPos1) && isempty(xPos2))

    % Increase thresold bounds if empty
    upT = upT + 0.5;
    lowT = lowT-0.5;

    % take the two indicies closes to the mid point
    Locat = find( mid < upT & mid > lowT);
    
    midX = floor(length(xVec)/2);
    
    % Lower Val
    xPos1 = xVec(max(Locat(Locat < midX)));
    xPos2 = xVec(min(Locat(Locat > midX)));
    idx = idx+1;

    if idx == 100
        error('Stuck in while loop in CR_getFWHM');
    end

end


% figure;
% plot(xVec, yVec)
% xline(xPos1,'-.') %
% xline(xPos2,'-.')


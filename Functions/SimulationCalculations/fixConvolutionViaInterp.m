function outputVec = fixConvolutionViaInterp(inputVec,interpMethod,num2interp)
% interpolate the first and last 5 values as they get cut off in the convolution.

xq = 1:length(inputVec);
outputVec = interp1(xq(num2interp:end-num2interp),inputVec(num2interp:end-num2interp),xq,interpMethod,'extrap');

% figure;
% plot(real(inputVec))
% hold on
% plot(real(outputVec))
% legend

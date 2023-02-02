function    [lowerLimit,highLimit] = updateFitLimits(cVal, lowerLimit, highLimit, idx)

% get an issue with some of these parameters changing non-linearly in
% realistic since. So try and move the limits more slowly if needed.

L1 = cVal(2)- 0.85*(cVal(2)-lowerLimit(idx));
L2 = lowerLimit(idx) + 0.1* lowerLimit(idx);
lowerLimit(idx) = min([L1,L2]);

H1 = cVal(2)+ 0.85*(highLimit(idx)-cVal(2));
H2 = highLimit(idx) - 0.1* highLimit(idx);
highLimit(idx) = max([H1,H2]);




% % use the fit values to update the fit.
% 
% % Taken from here( https://www.mathworks.com/matlabcentral/answers/499451-gaussian-fit-to-xy-data-and-extracting-fwhm) 
% % we can get the FWHM as:
% 
% FWHM = 2*sqrt(log(2))*cVal(3);
% 
% newLow = cVal(2) - 0.5*FWHM;
% newHigh = cVal(2) + 0.5*FWHM;
% 
% % Want to converge, so do not allow the limits to move further outwards
% if newLow < lowerLimit(idx) || newLow > highLimit(idx)
%     newLow = lowerLimit(idx);
% end
% if newHigh > highLimit(idx) || newHigh < lowerLimit(idx)
%     newHigh = highLimit(idx);
% end
% 
% lowerLimit(idx) = newLow;
% highLimit(idx) = newHigh;
function  std_resid = CR_calc_std_residuals( Xval2fit, Yval2fit, fitCoeffMat)

N_vals = length( Yval2fit ); % for normaliztion

Val_mean = mean( Yval2fit ,'omitnan');
Val_devsq = sum( (Yval2fit - Val_mean).^2); % sum of the squared difference from mean.
Val_lev = (1/N_vals) + ( (Yval2fit - Val_mean).^2 )./  Val_devsq; % should get a vector length (val2fit) for leverage

std_resid = zeros( length(fitCoeffMat),  1);

[x, y] = size(fitCoeffMat);

if min([x,y]) == 1
    y_fit = polyval(fitCoeffMat ,Xval2fit); % calculate fit curve
    resid = (y_fit-Yval2fit); % residuals between fit curve and data points
    s_Error = sqrt( (1/(N_vals-1)) * sum(resid.^2) ); % standard error
    std_resid = sum( resid ./ sqrt(s_Error .* (1- Val_lev) ) ); % standardized residuals
else
    for i = 1:length(fitCoeffMat) 
        
        y_fit = polyval(fitCoeffMat(i,:) ,Xval2fit); % calculate fit curve
        resid = (y_fit-Yval2fit); % residuals between fit curve and data points
        s_Error = sqrt( (1/(N_vals-1)) * sum(resid.^2) ); % standard error
        std_resid(i) = sum( resid ./ sqrt(s_Error .* (1- Val_lev) ) ); % standardized residuals
    end
end
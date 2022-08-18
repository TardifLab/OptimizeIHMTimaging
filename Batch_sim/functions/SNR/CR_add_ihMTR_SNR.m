function simResults = CR_add_ihMTR_SNR(simResults, noiseLvl)
% to be used with plotSimResultsIHMT_figures_SNR_v1.m

if nargin == 1   
  noiseLvl = 0.0005; 
end

% We just need a few columns from simResults:
% ihMT_proxy (10), and reference signal (17)
S0 = simResults(:,3);
ihMT_proxy = simResults(:,1) - simResults(:,2);

a= noiseLvl.^2 ./ (4 * S0.^2); % contribution to SNR from pos image
b= noiseLvl.^2 ./ (4 * S0.^2); % contribution to SNR from neg image
c= noiseLvl.^2 ./ (S0.^2);       % contribution to SNR from dual image
d= ((2*ihMT_proxy)./ (2 * S0.^2)).^2.*noiseLvl.^2 ; % contribution to SNR from S0 image
% Note in d, it should be negative ihMT proxy, but the square makes it not
% matter.

ihMTR = simResults(:,6);
ihMTR_SNR = ihMTR./sqrt( a+b+c+d);
ihMTR_SNR_eff = ihMTR_SNR ./simResults(:,10); % efficiency, where acquisition time is in column 10.

% Write out this way to prevent overwriting wrong thing

simResults(:,11) = ihMTR_SNR;
simResults(:,12) = ihMTR_SNR_eff;












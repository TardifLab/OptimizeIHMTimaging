function output = CR_add_MTR_SNR(simResults, noiseLvl)
% to be used with CNR_comp.m

%simResults=GM_simResults;

if nargin == 1   
  noiseLvl = 0.0005; 
end


% We just need a few columns from simResults:
% ihMT_proxy (10), and reference signal (17)
S0 = simResults(:,3);
sMTR_diff = S0 - simResults(:,1);

a= noiseLvl.^2 ./ (4 * S0.^2); % contribution to SNR from pos image
b= noiseLvl.^2 ./ (4 * S0.^2); % contribution to SNR from neg image
c= noiseLvl.^2 ./ (S0.^2);       % contribution to SNR from dual image
d= ((2*sMTR_diff)./ (2 * S0.^2)).^2.*noiseLvl.^2 ; % contribution to SNR from S0 image
% Note in d, it should be negative ihMT proxy, but the square makes it not
% matter.

sMTR = sMTR_diff./S0;
sMTR_SNR = sMTR./sqrt( a+b+c+d);



%%%%%%%%%%%% Repeat for the Dual MTR.

dMTR_diff = S0 - simResults(:,2);
d= ((2*dMTR_diff)./ (2 * S0.^2)).^2.*noiseLvl.^2 ; % contribution to SNR from S0 image
dMTR = dMTR_diff./S0;
dMTR_SNR = dMTR./sqrt( a+b+c+d);

output = [simResults sMTR_SNR dMTR_SNR];




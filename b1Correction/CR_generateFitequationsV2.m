function [fit_SS_eqn, fit_SS_eqn_sprintf, fit_SSsat, numTerms] = CR_generateFitequationsV2(M0b, b1, Raobs, MTsat_sim, threshold)

% V2 updates to let you easily set a flexible number of coefficients.
% use loops to write out equation. 

%figure; imshow3Dfull(MTsat_sim, [0 max(MTsat_sim,[],'all') ],jet)

if ~exist('threshold', 'var')
    threshold = 0.7;
end


B1_coeff_number = 5;
M0b_coeff_number = 2;
R1_coeff_number = 4;

%% Set up variables
[ M0b_mesh,b1_mesh,Raobs_mesh] = meshgrid(M0b',b1',Raobs'); % note: meshgrid swaps the first two dimensions, so we enter reverse order. 

MTsat_sim( MTsat_sim > threshold) = NaN;

M0b_l = M0b_mesh(:); % convert to 1x vector for fitting
b1_l = b1_mesh(:);
Raobs_l = Raobs_mesh(:);
fitz = MTsat_sim(:);

% Have issues if you have 0 MTsat values (from nan above), so remove
M0b_l( isnan(fitz)) = []; 
b1_l ( isnan(fitz)) = [];
Raobs_l ( isnan(fitz) ) = [];
fitz ( isnan(fitz)) = [];

%% Build the model and fit

% Generate the tricubic model. 
modelTerms = zeros( (1+B1_coeff_number)*(1+M0b_coeff_number)*(1+R1_coeff_number),3);
idx = 1;
for i = 0:M0b_coeff_number
    for j = 0:B1_coeff_number
        for k = 0:R1_coeff_number
            modelTerms(idx,1) = i;
            modelTerms(idx,2) = j;
            modelTerms(idx,3) = k;
            idx = idx+1;
        end
    end
end

numTerms = size(modelTerms,1);

tic
fit_SSsat = polyfitn([M0b_l, b1_l, Raobs_l], fitz, modelTerms); % ~ 5 seconds. 
toc

fit_SSsat.VarNames = {'M0b','b1','Raobs'}; % these will be the variable names required to evaluate them

%% write out equation
fit_SS_eqn = [];
fit_SS_eqn_sprintf = [];

idx = 1;
for i = 0:M0b_coeff_number
    for j = 0:B1_coeff_number
        for k = 0:R1_coeff_number
            
            fit_SS_eqn = strcat( fit_SS_eqn, ' + ',num2str(fit_SSsat.Coefficients(idx)),...
                '*M0b.^',num2str(i),'.*b1.^',num2str(j),'.*Raobs.^',num2str(k));
            
            fit_SS_eqn_sprintf = strcat( fit_SS_eqn_sprintf, ' + ',num2str(fit_SSsat.Coefficients(idx)),...
                '*M0b.^',num2str(i),'.*b1.^',num2str(j),'.*(%f).^',num2str(k));
            
            
            idx = idx+1;
        end
    end
end




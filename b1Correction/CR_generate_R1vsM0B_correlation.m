%% Function is meant to generate the plots for correlating M0B and R1obs
function [fitValues] = CR_generate_R1vsM0B_correlation( R1, M0b_map, mask, fitValues, outputImage, outputFitValues)

% fitValues = fitValues structure needing to be filled
% outputImage = filename for output correlation plot
% outputFitValues = filename for output fitValues structure

% Make sure R1 is in units 1/sec (not 1/ms!)

textClr = [0 0 0];

% Adjust mask in case M0b wasn't fit to all slices
mask(M0b_map <= 0) = 0;

% mask data
R1_p = R1(mask>0);
M0b_p = M0b_map(mask>0);

contrast_fit = zeros(1,2);

% Generate Fits
ft = fittype('poly1');
tmp = M0b_p;
tmp_r1 = R1_p; % convert from ms to sec
tmp_r1(tmp==0) = []; % to clean up the plots
tmp(tmp==0) = [];
M0b_fit = fit(tmp_r1,tmp,ft);
[R,P]= corrcoef([tmp, tmp_r1],'Alpha',0.05,'Rows','complete') ;
contrast_fit(1,1) = R(2,1);
contrast_fit(1,2) = P(2,1);
fitvals = coeffvalues(M0b_fit);
        
% Get plot limits
x_l = min(R1_p); 
x_h = max(R1_p); 
y_h = 0.25; %max(M0b_p);  % Set manually...

% Generate Plots       
f1 = figure;
heatscatter(tmp_r1,tmp, spring)
xlim([x_l x_h])
ylim([0 y_h])
hold on
hline = refline(fitvals(1),fitvals(2));     hline.Color = 'blue';     hline.LineWidth= 2;
caption = sprintf('M_{0,app}^B = %.2g * R_{1,obs} %.2g', fitvals(1), fitvals(2));
text(1.15*x_l, 0.9*y_h, caption, 'FontSize', 16, 'Color', textClr,'fontweight', 'bold');    
caption2 = sprintf('r = %.2f', contrast_fit(1,1));
text(1.15*x_l, 0.8*y_h, caption2, 'FontSize', 16, 'Color', textClr,'fontweight', 'bold');
ax = gca;
ax.FontSize = 20; 
xlabel('R_{1,obs} (1/s)', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('M_{0,app}^B', 'FontSize', 20, 'FontWeight', 'bold')
colorbar('off')
legend('hide')
saveas( f1, outputImage ) 
    
  
% Export the results to fitvalues and save
fitValues.Est_M0b_from_R1 = strcat( num2str(fitvals(1) ),' *Raobs + ', num2str(fitvals(2) ) );
save(outputFitValues,'fitValues')


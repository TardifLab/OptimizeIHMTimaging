%% Function is meant to generate the plots for correlating M0B and R1obs
function [fitValues] = CR_generate_R1vsM0B_correlation( R1, M0b_map, mask,...
                    fitValues, outputImage, outputFitValues, fitDegree)

% fitValues = fitValues structure needing to be filled
% outputImage = filename for output correlation plot
% outputFitValues = filename for output fitValues structure

% Make sure R1 is in units 1/sec (not 1/ms!)

% New - fitDegree:
% Specify 1 or 2 as the polynomial order to fit.

% Written by Christopher Rowley 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist( "fitDegree", "var")
    fitDegree = 1;
end

textClr = [0 0 0];

% Adjust mask in case M0b wasn't fit to all slices
mask(M0b_map <= 0) = 0;

% mask data
R1_p = R1(mask>0);
M0b_p = M0b_map(mask>0);

% Get plot limits
x_l = min(R1_p); 
x_h = max(R1_p); 
y_h = 0.25; %max(M0b_p);  % Set manually...

% Prep data:
tmp = M0b_p;
tmp_r1 = R1_p; % convert from ms to sec
tmp_r1(tmp==0) = []; % to clean up the plots
tmp(tmp==0) = [];

% Start the plot: 
 
f1 = figure;
heatscatter(tmp_r1,tmp, spring)
xlim([x_l x_h]); 
ylim([0 y_h]); 
hold on

if fitDegree == 1
    
    % Generate Fit
    [M0b_fit, gof] = fit( tmp_r1, tmp, 'poly1');
    fitvals = coeffvalues(M0b_fit);

    % Generate plot text    
    caption = sprintf('M_{0,app}^B = %.2g * R_{1,obs} %.2g', fitvals(1), fitvals(2));
    caption2 = strcat( 'R^2 = ', sprintf('%.2f', gof.rsquare ) );
    
    % Export the results to fitvalues
    fitValues.Est_M0b_from_R1 = strcat( num2str(fitvals(1) ),...
        ' *Raobs + ', num2str(fitvals(2) ) );

elseif fitDegree == 2
    
    % Linear 'should' work. This might come up with other issues with
    % spoiling or inaccurate B1 mapping
        
    [M0b_fit, gof] = fit( tmp_r1, tmp, 'poly2');
    fitvals = coeffvalues(M0b_fit);

    caption = sprintf('M_{0,app}^B = %.2g * R_{1,obs}^2  %.2g * R_{1,obs}+ %.2g',...
                                fitvals(1), fitvals(2), fitvals(3));
    caption2 = strcat( 'R^2 = ', sprintf('%.2f', gof.rsquare ) );
    
    % Export the results to fitvalues and save
    fitValues.Est_M0b_from_R1 = strcat( num2str(fitvals(1) ),...
        ' *Raobs.^2 + ',  num2str(fitvals(2) ),' *Raobs +', num2str(fitvals(3) ) );
    
else
    error('Polynomials of degree 1 and 2 are only supported at this time. If this doesnt fit, check your simulations');
end

% Finish plot with common code:
hline = plot( M0b_fit );  hline.Color = 'blue';     hline.LineWidth= 2;
text(1.05*x_l, 0.9*y_h, caption, 'FontSize', 14, 'Color', textClr );    
text(1.05*x_l, 0.8*y_h, caption2, 'FontSize', 14, 'Color', textClr);
ax = gca;    ax.FontSize = 18; 
xlabel('R_{1,obs} (1/s)', 'FontSize', 18, 'FontWeight', 'bold')
ylabel('M_{0,app}^B', 'FontSize', 18, 'FontWeight', 'bold')
colorbar('off')
legend('hide')
saveas( f1, outputImage ) 
    

save(outputFitValues,'fitValues')

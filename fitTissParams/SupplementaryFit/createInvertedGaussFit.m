function cVal = createInvertedGaussFit(xData, yData)

% Generated using CFtool, then modified.
% set X guess as starting guess for FWHM. gets around the different orders
% of magnitude we face with parameters.

% Need to flip data to get a maximum instead of minimum:
mxY = max(yData);
yData = mxY + (mxY - yData);

setA = mxY*2; % peak should have been zero, fliped up.

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.Lower = [setA xData(1) 0];
% opts.StartPoint = [setA xData(2) xData(2)]; 
% opts.Upper = [setA xData(3) Inf];

opts.Lower = [setA*0.95 xData(1) 0];
opts.StartPoint = [setA xData(2) xData(2)]; 
opts.Upper = [setA*1.05 xData(3) Inf];


% Fit model to data.
[fitresult, ~] = fit( xData, yData, ft, opts );
cVal = coeffvalues(fitresult);

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'y vs. x', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'x', 'Interpreter', 'none' );
% ylabel( 'y', 'Interpreter', 'none' );
% grid on



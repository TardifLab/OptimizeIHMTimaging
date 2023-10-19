function CR_generateFitSimFiguresRelB1(M0b, b1, RaobsV, MTsat_sim, fit_SS_eqn, fileName)

% Modified version where the B1 axis is just the relative B1+ field. 
% To be used with newer versions using the saturation flip angle instead of
% B1rms. 

% Note the values M0b, b1, and Raobs need to be labelled as such for the
% proper evaluation of fit_SS_eqn. Thus RaobsV, is the Vector

% Also modified the 'slice' to find the index closest to 1000ms;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ M0b_mesh,b1_mesh,Raobs_mesh] = meshgrid(M0b',b1',RaobsV'); % note: meshgrid swaps the first two dimensions, so we enter reverse order. 

M0b =  M0b_mesh;
b1 = b1_mesh;
Raobs = Raobs_mesh;
z = eval(fit_SS_eqn);


%% Generate figure to show result at T1 = 1000ms;
[ ~, slice ] = min( abs( RaobsV - 1) );

slice_b1 = squeeze(b1_mesh(:,:,slice));
slice_M0b = squeeze(M0b_mesh(:,:,slice));
slice_MT = squeeze(MTsat_sim(:,:,slice));
slice_z = squeeze(z(:,:,slice));


f1 = figure;
surf(slice_b1, slice_M0b, slice_MT,'FaceAlpha',0.5)
hold on
surf(slice_b1, slice_M0b, slice_z,'FaceAlpha',0.5)
ax = gca; ax.FontSize = 16; 
xlabel('Relative B_1^+', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('M_{0b}', 'FontSize', 16, 'FontWeight', 'bold')
zlabel('MT_{sat}', 'FontSize', 16, 'FontWeight', 'bold')
%title('8k Dual')
%legend('sim','meas')

saveas(f1,fileName)

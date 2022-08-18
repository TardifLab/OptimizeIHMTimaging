function [Beta1, Beta2] = ComputeDiffusionBetas(M_in, Params, gam, G, G_t)
% As in:
% Jochimsen, T.H., Schäfer, A., Bammer, R., Moseley, M.E., 2006. Efficient 
% simulation of magnetic resonance imaging with Bloch-Torrey equations using 
% intra-voxel magnetization gradients. J. Magn. Reson. 180, 29–38. https://doi.org/10.1016/j.jmr.2006.01.001
 
% Slight change to the equations used to pull out time and insert 'D'

% gam = 2*pi*42.577478518e6; % gyromagnetic ratio in rad*T/s 

N_spin = size(M_in,2);
if G > 0
    % rad/T/s  * T/m * delta_m /m = % delta radians per second per meter
    GradOmega = gam*G*( Params.ReadoutResolution/ N_spin)/( Params.ReadoutResolution/ N_spin); % change in radians per meter
    h = G_t*GradOmega; % rad/s per meter* s = rad per meter
    gradPhi =  ComputeGradPhi(M_in);
else
    h = 0;
    gradPhi = 0;
end



% Beta 1 is longitudinal component, only spin diffusion
% Beta 2 is transverse, with potential gradient spoiling added

% Diffusion Terms - combine xyz
% [d] = ComputeIsoDiffusionTerm(M_in, Params);
% Beta1 = d*Params.D;
% Beta2 = Params.D*( 1/3*h.^2 - h.*gradPhi + d);


% Diffusion Terms - separate transverse and longitudinal
[d_perp, d_z] = ComputeIsoDiffusionTerm(M_in, Params);
Beta1 = d_z*Params.D; 
Beta2 = Params.D*( 1/3*h.^2 - h.*gradPhi + d_perp);


%% We get instability here, and considering the formalism presented in Jochimsen,
% we should be able to just take the median value of Beta 2, and Beta1 is
% just a single number anyway.

% Beta1 = median(Beta1);
% Beta2 = median(Beta2);





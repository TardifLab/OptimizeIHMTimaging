function M_final = ApplySpinEvolution_v2( Params, M_in, t, MTC, Gradient )

% Relaxation, gradient spoiling and diffusion effects with bound pool exchange 
% Diffusion limited to sections longer than 10ms due to instability below
% this

% Version 2 converts transverse magnetization to a complex number for
% computing precession from gradient, and diffusion.

% Input is:
% Params structure that stores a bunch of variables
% M_t is a 5xnSpins vector
% t is the time over which this step occurs

% Need values in Params structure for:
% Params.G_t -> gradient time in seconds
% Params.G -> gradient height in mT/m.
% Params.D -> diffusion coefficient in m/s!


% Relevant sources:
% Gudbjartsson, H., Patz, S., 1995. NMR Diffusion Simulation on Conditional Random Walk. IEEE Trans. Med. Imaging 14, 636–642.
% Kiselev, V.G., 2003. Calculation of diffusion effect for arbitrary pulse sequences. J. Magn. Reson. 164, 205–211. https://doi.org/10.1016/S1090-7807(03)00241-6
% Yarnykh, V.L., 2010. Optimal radiofrequency and gradient spoiling for improved accuracy of T1 and B1 measurements using fast steady-state techniques. Magn. Reson. Med. 63, 1610–1626. https://doi.org/10.1002/mrm.22394
% Jochimsen, T.H., Schäfer, A., Bammer, R., Moseley, M.E., 2006. Efficient simulation of magnetic resonance imaging with Bloch-Torrey equations using intra-voxel magnetization gradients. J. Magn. Reson. 180, 29–38. https://doi.org/10.1016/j.jmr.2006.01.001
% Sled, J.G., Pike, G.B., 2000. Quantitative Interpretation of Magnetization Transfer in Spoiled Gradient Echo MRI Sequences. J. Magn. Reson. 145, 24–36. https://doi.org/10.1006/jmre.2000.2059
% Kose, R., Kose, K., 2017. BlochSolver: A GPU-optimized fast 3D MRI simulator for experimentally compatible pulse sequences. J. Magn. Reson. 281, 51–65. https://doi.org/10.1016/j.jmr.2017.05.007

I = eye(5);
N_spin = size(M_in,2);
gam = 2*pi*42.577478518e6; % gyromagnetic ratio in rad*T/s 


M_perp = complex( M_in(1,:), M_in(2,:));
M_z = M_in(3,:);

%% Toggle between two spoiling, one for MT and one for water:
if Gradient
    
    if MTC % Spoil Sat pulse
        G_t = Params.G_t_MT;
        G = Params.GradientSpoilingStrength_MT/1000;
    
    else % Spoil water excitation
        G_t = Params.G_t;
        G = Params.GradientSpoilingStrength/1000;
    end
    
    %% Do gradient dephasing first:
    dis = linspace(0,Params.ReadoutResolution, N_spin); % meters
    % theta = gam* G_t* G*dis ; % seconds * rad/T * T/m *m -> radians
    
    M_perp = M_perp.*exp(1i*gam*G*dis*G_t); 

else % if you only want diffusion
    G_t = 0;
    G = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Diffusion 
% Added a restriction to not compute over very short time scales.
if Params.ModelSpinDiffusion && (t >= 3e-3)
    [M_perp2, M_z2] = ComputeSpinDiffusion(M_perp, M_z, Params, gam, G, G_t, t);
else
    M_perp2 = M_perp;
    M_z2 = M_z;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% relaxation terms.
M_out = M_in;
M_out(1,:) = real(M_perp2);
M_out(2,:) = imag(M_perp2);
M_out(3,:) = M_z2;

E = [-1/Params.T2a, 0,  0,     0,      0;...
      0, -1/Params.T2a,  0,     0,      0;...
      0, 0, -Params.kf-Params.Ra, Params.kr,   0;...
      0, 0, Params.kf, -Params.kr-Params.R1b,   0;...
      0, 0,         0,    0, -1/Params.T1D];

B = [0 0 Params.Ra*Params.M0a, Params.R1b*Params.M0b, 0]';

M_final = zeros(size(M_out));

for i = 1:N_spin
        M_final(:,i) =expm(E*t)*M_out(:,i) + (expm(E*t) - I)* (E\B);
end


% figure
% plot(sqrt(sum(M_in(1:2,:).^2)))
% hold on; 
% plot(sqrt(sum(M_final(1:2,:).^2)))
% legend
% 
% figure
% plot((M_in(2,:)))
% hold on; 
% plot((M_out(2,:)))
% plot((M_in(1,:)))
% hold on; 
% plot((M_out(1,:)))
% legend

% figure
% plot(M_in(3,:))
% hold on; 
% plot(M_final(3,:))
% legend














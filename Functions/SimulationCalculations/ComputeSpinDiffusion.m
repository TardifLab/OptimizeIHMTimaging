function [M_perp2, M_z2] = ComputeSpinDiffusion(M_perp, M_z, Params, gam, G, G_t, t)

% Based on convolution method from:
% Gudbjartsson, H., Patz, S., 1995. NMR Diffusion Simulation on Conditional Random Walk. IEEE Trans. Med. Imaging 14, 636–642.

% According to their framework, they:
% 1. Precess spins as if static; [Done prior to this function]

% ********************************************************************
% 2. Convolve the magnetization with the correction kernel c(x,t) to
% account for diffusion
% ********************************************************************


% gam = 2*pi*42.577478518e6; % gyromagnetic ratio in rad*T/s 
N_spin = size(M_perp,2);


% We will assume gradient is constant over the duration. If not, then
% modify G to get equival G strength.

if G_t ~= t
    % modify G
    G = G*G_t/t;
end

x = linspace(-Params.ReadoutResolution/2,Params.ReadoutResolution/2, N_spin); % meters

c1 = 1/sqrt(2*Params.D*t*2*pi);
c2 = -1*x.^2/(4*Params.D*t);
c3 = -gam^2*G^2*Params.D*t^3/12;
c4 = 1i*gam*G*x./2*t;

c = c1*exp(c2+c3+c4);

% Normalize for Prob Dist Function
c = c./trapz(c);

M_perp2 = conv(M_perp, c,'same');
M_perp2 = fixConvolutionViaInterp(M_perp2,'pchip',6);

% Repeat for longitudinal magnetization:
c = c1*exp(c2);
c = c./trapz(c);

M_z2 = conv(M_z, c,'same');
M_z2 = fixConvolutionViaInterp(M_z2,'linear',10);


% figure;
% plot(real(M_perp))
% hold on
% plot(real(M_perp2))
% legend
% 
% figure;
% plot(real(M_z))
% hold on
% plot(real(M_z2))
% legend
































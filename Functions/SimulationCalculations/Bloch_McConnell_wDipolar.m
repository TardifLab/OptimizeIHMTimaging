function E_rf = Bloch_McConnell_wDipolar( Params, delta, w1)
%% Requires functions tied to qMRlab for computing lineshapes


% Keep track of Water X,Y,Z as the simulations use transverse for signal
% calculation. It will be important for the direct effect. 

% Solving equations as written out in:
% Henkelman, R.M., Huang, X., Xiang, Q. ‐S, Stanisz, G.J., Swanson, S.D., Bronskill, M.J., 1993.
%   Quantitative interpretation of magnetization transfer. Magn. Reson. Med. 29, 759–766. https://doi.org/10.1002/mrm.1910290607

% Sled JG, Pike GB. Quantitative Interpretation of Magnetization Transfer in 
%   Spoiled Gradient Echo MRI Sequences. J Magn Reson. 2000;145(1):24-36. doi:10.1006/jmre.2000.2059

% With the added dipolar order from: 
% Morrison C, Stanisz G, Henkelman RM. Modeling Magnetization Transfer for 
%   Biological-like Systems Using a Semi-solid Pool with a Super-Lorentzian
%   Lineshape and Dipolar Reservoir. J Magn Reson Ser B. 1995;108(2):103-113. 
%   doi:10.1006/jmrb.1995.1111

% Modification to dual alternate and single saturation values follows from:
% Lee, J.S., Khitrin, A.K., Regatte, R.R., Jerschow, A., 2011. Uniform 
%   saturation of a strongly coupled spin system by two-frequency irradiation. 
%   J. Chem. Phys. 134, 1–6. https://doi.org/10.1063/1.3600758

% This is perhaps made more clear in:
% Manning, A.P., Chang, K.L., MacKay, A.L., Michal, C.A., 2017. The physical
%   mechanism of "inhomogeneous" magnetization transfer MRI. J. Magn. Reson. 
%   274, 125–136. https://doi.org/10.1016/j.jmr.2016.11.013

%% Calculate some other parameters:
% kf = (Params.R*Params.M0b);
% kr = (Params.R*Params.M0a);
R2a = 1/Params.T2a;

%w1 = 2*pi *42.577478518 * b1;  % assume B1 is in microTesla, and drop the 10^6 from gamma. w1 in rad/s

%% For consistency, use computeG code from qMRLab
% however, this requires us to use abs(delta)
G = computeG( abs(delta), Params.T2b, Params.lineshape);

if strcmp(Params.lineshape, 'SuperLorentzian') % default is gaussian
    Rrfb = pi.*w1.^2.* G;
    wloc = sqrt(1/(15*Params.T2b^2)); % saturation impacting dipolar pool (Morrison et al 1995)
    
elseif strcmp(Params.lineshape, 'Lorentzian') % default is gaussian
    Rrfb = pi.*w1.^2.* G;
    wloc = sqrt(1/(3*Params.T2b^2)); % saturation impacting dipolar pool (Morrison et al 1995)
    
else % default is 'Gaussian'
    Rrfb = w1.^2 .* G;
    wloc = sqrt(1/(3*Params.T2b^2)); % saturation impacting dipolar pool (Morrison et al 1995)
end


% These need to be calculated after the lineshape details are sorted
Omega = 2*pi*delta/wloc;

if ~Params.IncludeDipolar
    % hack to have super fast relaxation, and will set omega to 0
    % permit quick check against two pool model results.
    Params.T1D = 1e-6; 
    Omega = 1;
end


kf = Params.kf;
kr = Params.kr;

if delta > 0 && strcmp( Params.freqPattern,'dualContinuous')
    % Provotorov theory, differential saturation for simultaneous dual
    E_rf =[ -R2a,   -2*pi*delta,    0,       0,       0; ...       % Water X
            2*pi*delta,     R2a,  -w1,       0,       0;...        % Water Y
             0,          w1, -Params.Ra-kf,  kr,      0;...        % Water Z
             0,           0,       kf, -Rrfb-kr-Params.R1b,   0; ...% Bound Z
             0,           0,        0,       0, -1/Params.T1D ];   % Dipolar

else % Water excitation
    E_rf =[ -R2a,  -2*pi*delta,     0,           0,           0; ...          % Water X
          2*pi*delta,  -R2a,      -w1,           0,           0;...           % Water Y
             0,          w1, -Params.Ra-kf,     kr,           0;...           % Water Z
             0,           0,       kf, -Rrfb-kr-Params.R1b,   Omega*Rrfb; ...  % Bound Z
             0,           0,        0,  Omega*Rrfb, -(Omega^2)*Rrfb-1/Params.T1D ];% Dipolar
end
























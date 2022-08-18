function M_out = XYmag_Spoil(Params, M_in, t, MTC, Gradient)

% Intended to spoil the magnetization at the end of a magnetizaton
% saturation pulse train.

% 'Gradient' is a binary variable


% M_in is a 5xN matrix, where N is the number of spins. 

if Params.PerfectSpoiling
    M_out = SpinEvolution_Relaxation( Params, M_in, t); % Calculate relaxation over period t
    M_out(1:2,:) = 0; % Spoil. Easier to do if you are keeping isochromats
else

    M_out = ApplySpinEvolution_v2( Params, M_in, t, MTC, Gradient );

end

% can get issues here if there is no difference.
M_out(isnan(M_out)) = M_in(isnan(M_out));
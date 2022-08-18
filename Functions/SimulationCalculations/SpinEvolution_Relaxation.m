function M_out = SpinEvolution_Relaxation( Params, M_in, t)

% Calculate spin evolution in the absence of RF and gradients

% Input is:
% Params structure that stores a bunch of variables
% M_in is a 5x Number_spin vector
% t is the time over which this step occurs

% You have the option to insert a wide matrix with lots of spins, or pool
% all the spins, and calculate 1 vector. 

ns = size(M_in,2);

% kf = (Params.R*Params.M0b);
% kr = (Params.R*Params.M0a);
B = [0 0 Params.Ra*Params.M0a, Params.R1b*Params.M0b, 0]';
I = eye(5);

% Evolution Matrix
E = [-1/Params.T2a, 0,  0,     0,      0;...
      0, -1/Params.T2a,  0,     0,      0;...
      0, 0,        -Params.kf-Params.Ra, Params.kr,   0;...
      0, 0,        Params.kf, -Params.kr-Params.R1b,   0;...
      0, 0,         0,    0, -1/Params.T1D];


M_out = zeros(5,ns);
for i = 1:ns
    M_out(:,i) = expm(E*t) * M_in(:,i) + (expm(E*t) - I)* (E\B);
end



% M_in = [0.3; 0.3; 0.8; 0.02; 0.1];
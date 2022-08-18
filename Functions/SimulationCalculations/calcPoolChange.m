function [M_out, t_out] = calcPoolChange(A_mat, B_mat, I, t, M_in, t_in)

% This is intended to clean up code to calculate (n) pool changes. 
% Inputs:
% Precalculated A matrix containing RF and relaxation info (size nxn)
% Preset B relaxation matrix (size nx1)
% I = eye(length(B)) -> calculate outside for performance
% t = time duration of the relaxation step
% M_1 = the input magnetization to this time step
% t_1 = the input time - only used to update time vector for plotting

% Output
% M_2 = output magnetization
% t_2 = new sequence time - used for plotting

% Built to handle multiple isochromats if needed
% ns = size(M_in,3);
% M_out = zeros(size(M_in));
% AExp = expm(A_mat*t);
% for z = 1:ns
%     M_out(:,z) = AExp * M_in + (AExp - I)* (A_mat\B_mat); % Update Magnetization.
% end
% t_out = t_in + t;


AExp = expm(A_mat*t);
AEnd = (AExp - I)* (A_mat\B_mat);
M_out = pagemtimes(AExp, M_in) + AEnd;

t_out = t_in + t;

function [Output_phase, last_inc] = IncrementRFspoilPhase( input_phase, Params, last_inc )

% Use the quadratic phase increment scheme. 
% In general, value is based on vendor GE = 117, Siemens = 50, Philips = 150 
% 𝜙n = 𝜙0 /2 * (n + 1) n 

%% For Siemens, precompute for the whole excitation train:
inc = Params.RFspoilingPhaseInc; % increment in degs

Output_phase = zeros(Params.numExcitation,1);
% Output_phase(1) = input_phase;

for i = 1:Params.numExcitation
    last_inc = inc + last_inc;
    if i == 1
        Output_phase(i) = input_phase + last_inc;
    else
        Output_phase(i) = Output_phase(i-1) + last_inc;
    end
end

% take the remainder so it is within 360 deg
Output_phase = rem(Output_phase, 360);

% Convert to rad: -> nevermind, easiest to keep track in degrees
% Output_phase = pi*Output_phase./ 180;

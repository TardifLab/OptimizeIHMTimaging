%% This code is meant to solve the signal equations presented in Deichmann and Haase for FLASH imaging
% This requires knowledge of:
% - flip angle (in degrees)
% - repetition time TR

% Also need estimates of:
% B1field - in relative units, where 1 == nominal flip angle
% M0 - apparent signal
% T1- longitudinal relaxation rate


function Sig= CR_FLASH_solver(flipAngle, TR, B1field, M0, T1)

flip_a = (flipAngle*B1field) * pi / 180; % correct for B1 and convert to radians
x = cos(flip_a) ;
y = exp(-TR/T1);

% Solve for magnetization
M = M0*(1-y) ./ (1-x*y);

% We read out the value M with sin(flip)
Sig = sin(flip_a) * M;




% % DERIVE M = M0*(1-y) ./ (1-x*y);
% M = magnetization before flip angle; 
% M1 = magnetization after flip
% M2 = magnetization after T1 recovery
% @ steady state M2 == M

% M2 = M0(1-exp(-t/T1)) + M1*exp(-t/T1)
% M1 = M*cos(a)
% Substitute
% M2 = M = M0(1-exp(-t/T1)) + (M*cos(a))*exp(-t/T1)
% Convert to the x and y values above for easier reading
% M = M0(1-y) + (M*x*y)
% M - M*x*y = M0(1-y)
% M(1-x*y) = M0(1-y)
% M = M0(1-y) / (1-x*y)


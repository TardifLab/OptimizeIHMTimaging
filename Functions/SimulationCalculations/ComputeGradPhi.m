function gradPhi =  ComputeGradPhi(M_in)

% Compute the gradient of the spin phase

Mx = M_in(1,:);
My = M_in(2,:);
M_perp = Mx.^2 + My.^2;


% gradPhi = 0 where M_perp == 0;

% For M_perp ~= 0 ...
% v1 = Mx.* gradient(My, Params.ReadoutResolution./size(M_in,2));
% v2 = My.* gradient(Mx, Params.ReadoutResolution./size(M_in,2));

v1 = Mx.* gradient(My);
v2 = My.* gradient(Mx);

gradPhi = (v1 - v2) ./ M_perp;


% threshold for rounding issues
gradPhi = median(gradPhi);
% Note M_perp should be squared, which cancels the square root.

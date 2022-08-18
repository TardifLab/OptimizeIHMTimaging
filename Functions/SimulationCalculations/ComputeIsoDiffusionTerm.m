%function [d] = ComputeIsoDiffusionTerm(M_in, Params)

function [d_perp, d_z] = ComputeIsoDiffusionTerm(M_in, Params)


% as in Jochimsen, T.H., Schäfer, A., Bammer, R., Moseley, M.E., 2006.

% Based on their formalism, each isochromat is restrained by some sphere
% with radius |M_i|

Mi = sum(M_in(1:3,:).^2);
gradVal = gradient(M_in, Params.ReadoutResolution./size(M_in,2) ).^2;
%d = sum(gradVal(1:3,:),1)./ Mi ;
d_perp = sum(gradVal(1:2,:),1)./ Mi ;
d_z = gradVal(3)./ Mi ;


% May need to incorporate the measurement spacing for the gradient...
% if so, this gradient(M_in) becomes
% gradient(M_in, Params.ReadoutResolution./size(M_in,2) )

% Since it is multiplied by a diffusion coefficient term that has units m^2
% Here is the equivalent fully written out:
% Mi = sqrt(sum(M_in(1:3,:).^2));
% gradVal = gradient(M_in, Params.ReadoutResolution./size(M_in,2) ).^2;
% d = sum(gradVal,1)./ (Mi.^2) ;
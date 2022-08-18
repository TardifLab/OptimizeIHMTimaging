function R = RotationMatrix_withBoundPool_MP2RAGE(fa, ph, Params, invNumber)
% Rotation matrix for a flip angle and phase:
% Note that these matlab functions assume both values are in radians

if invNumber == 1
    Erfb = exp(-Params.Rrfb_exc1*Params.WExcDur);
elseif invNumber == 2
    Erfb = exp(-Params.Rrfb_exc2*Params.WExcDur);
else
    error('specify invNumber as 1 or 2')
end

Erfd = exp(-Params.Rrfd_exc*Params.InvPulseDur); % doing to be exp(0)

R = [cos(fa)+(1-cos(fa))*cos(ph)^2, (1-cos(fa))*sin(ph)*cos(ph),    -sin(fa)*sin(ph), 0, 0;...
    (1-cos(fa))*sin(ph)*cos(ph),    cos(fa)+(1-cos(fa))*sin(ph)^2 ,  sin(fa)*cos(ph), 0, 0;...
    sin(fa)*sin(ph),                     -sin(fa)*cos(ph),                   cos(fa), 0, 0;...
                  0,                                  0,                       0,  Erfb, 0;...
                  0,                                  0,                        0, 0, Erfd];















































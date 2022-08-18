function R = RotationMatrix_withBoundPool(fa, ph, Params)
% Rotation matrix for a flip angle and phase:
% Note that these matlab functions assume both values are in radians


Erfb = exp(-Params.Rrfb_exc*Params.WExcDur);
Erfd = exp(-Params.Rrfd_exc*Params.WExcDur);

R = [cos(fa)+(1-cos(fa))*cos(ph)^2, (1-cos(fa))*sin(ph)*cos(ph),    -sin(fa)*sin(ph), 0, 0;...
    (1-cos(fa))*sin(ph)*cos(ph),    cos(fa)+(1-cos(fa))*sin(ph)^2 ,  sin(fa)*cos(ph), 0, 0;...
    sin(fa)*sin(ph),                     -sin(fa)*cos(ph),                   cos(fa), 0, 0;...
                  0,                                  0,                       0,  Erfb, 0;...
                  0,                                  0,                        0, 0, Erfd];















































function PerpMag = TransverseMagnetizationMagnitude(M)
% calculate perpendicular magnetization from a group of isochromats

avgx = mean(M(1,:));
avgy = mean(M(2,:));
PerpMag = sqrt(avgx^2 +avgy^2);
    
%
function [cw, ccw, still] = headTurnGroups(angularHeadVelocity, p)
    v = abs(angularHeadVelocity);
    eliminateIdx = v > p.hdStillSpeedThreshold & v < p.hdTurnSpeedThreshold;

    angularHeadVelocity(eliminateIdx) = [];
    v(eliminateIdx) = [];

    cw = angularHeadVelocity(angularHeadVelocity > p.hdTurnSpeedThreshold);
    ccw = angularHeadVelocity(angularHeadVelocity < - p.hdTurnSpeedThreshold);
    still = angularHeadVelocity(v < p.hdStillSpeedThreshold);
end

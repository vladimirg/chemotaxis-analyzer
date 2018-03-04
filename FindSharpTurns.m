function sharp_turn_mask = FindSharpTurns(track, fps)
% Output is relative to velocity vector origins.
% So, if there are N steps, the output length will be N-1 steps, with the
% first value always being zero (as the turn state is undefined for the
% first velocity vector).
% Points at which a sharp turn are detected will be marked by
% SegmentMasks.SharpTurn, otherwise they will be marked by a 0.

% TODO: verify this, as some of our sharp turns are not sharp at all!
sharp_turn_angle = 50;

turning_rate = GetTurningRate(track, fps);
sharp_turn_mask = [0; abs(turning_rate(1:end-1)) > sharp_turn_angle;] * SegmentMasks.SharpTurn;

end
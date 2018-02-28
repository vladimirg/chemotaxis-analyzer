function sharp_turn_mask = FindSharpTurns(track, fps)
% Output is relative to velocity vector origins.
% So, if there are N steps, the output length will be N-1 steps, with the
% first value always being zero (as the turn state is undefined for the
% first velocity vector).

% TODO: verify this, as some of our sharp turns are not sharp at all!
sharp_turn_angle = 50;

turning_rate = GetTurningRate(track, fps);
sharp_turn_mask = [0; abs(turning_rate(1:end-1)) > sharp_turn_angle;] * SegmentMasks.SharpTurn;

end






function ixes = FindSharpTurnsOld(track, fps)
% Output is relative to velocity vector origins.
% So, if there are N steps, the output length will be N-1 steps, with the
% first value always being zero (as the turn state is undefined for the
% first velocity vector).

sharp_turn_angle = 50;

as = track.anglesBetweenSteps;
as_len = length(as);

% Kinks are defined as a single step in a sharp angle, but then the next
% step returning the original orientation. Thus, for every candidate, we
% check if it reverses the previous turn.

% TODO: a kink is made up of two steps, and theoretically there can be a
% series of kinks made up of an even number of individual kinks. This can
% be tricky to detect if the first/last step is a true sharp turn, and we
% are not handling this even/odd distinction at present.

candidate_ixes = (abs(as) > sharp_turn_angle) * SegmentMasks.SharpTurn;
[candidate_segments, seg_num] = SegmentVector(candidate_ixes);

% We only allow kinks as pairs flanked by normal angles:
for seg_ix = 1:seg_num
    if candidate_segments(seg_ix, 3) ~= SegmentMasks.SharpTurn
        continue
    end
    
    first_kink_ix = candidate_segments(seg_ix, 1);
    second_kink_ix = candidate_segments(seg_ix, 2);
    
    if mod(second_kink_ix - first_kink_ix, 2) == 1 % Even number of kinks
        angle_diff = mod(abs(sum(as(first_kink_ix:second_kink_ix))), 360);
        if  min(angle_diff, 360 - angle_diff) <= sharp_turn_angle
            candidate_ixes(first_kink_ix:second_kink_ix) = SegmentMasks.Kink;
        end
    else % Odd number of kinks
        % Look ahead
        look_ahead_kink = true;
        if second_kink_ix < as_len
            look_ahead_kink = false;
            look_ahead = mod(abs(sum(as(first_kink_ix:second_kink_ix+1))), 360);
            if min(look_ahead, 360 - look_ahead) <= sharp_turn_angle
                look_ahead_kink = true;
            end
        end
        
        % Look back
        look_back_kink = true;
        if first_kink_ix > 1
            look_back_kink = false;
            look_back = mod(abs(sum(as(first_kink_ix-1:second_kink_ix))), 360);
            if min(look_back, 360 - look_back) <= sharp_turn_angle
                look_back_kink = true;
            end
        end
        
        if look_ahead_kink || look_back_kink
            candidate_ixes(first_kink_ix:second_kink_ix) = SegmentMasks.Kink;
        end
    end
end

ixes = [0; candidate_ixes];

end
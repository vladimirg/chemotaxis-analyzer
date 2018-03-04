function fixed_bearings = CircularizeBearings(bearings)
% Bearings assume a value in the [-180, +180] interval. However, for some
% purposes (such as tracking the turning rate and turning bias) this is
% inconvenient, since bearing at t=N can be -170, bearing at t=N+1 can be
% +170, and the most likely explanation is that the worm turned -20 degrees
% (20 degrees clockwise). This function takes a list of bearings and
% "circularizes" them - in the above example, the bearing at t=N+1 will be
% changed to -190 degrees.

fixed_bearings = bearings;
added_degrees = 0;

for ix = 2:length(bearings)
    prev_ix = ix - 1;
    if bearings(prev_ix) > 90 && bearings(ix) < -90
        added_degrees = added_degrees + 360;
    elseif bearings(prev_ix) < -90 && bearings(ix) > 90
        added_degrees = added_degrees - 360;
    end
    fixed_bearings(ix) = bearings(ix) + added_degrees;
end

end
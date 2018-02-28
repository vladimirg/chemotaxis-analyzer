function fixed_bearings = FixBearings(bearings)

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
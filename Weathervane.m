function weathervane = Weathervane(track, movie_features)
% Compute the weathervane stats of a track. Per Iino, this is the change in
% bearing between the bearing of a local velocity vector and the bearing to
% the point 1 mm "into the future".
% The output is an Nx4 matrix, where the columns are:
% 1 - the bearing of the local velocity vector.
% 2 - the change in bearing, as defined above.
% 3 - the frame of the velocity vector.
% 4 - the distance between the origin of the velocity vector and the
%     center of the odorant.

segmented_track = SegmentTrackPath(track, movie_features.fps);
long_run = bitand(segmented_track, SegmentMasks.LongRun) == SegmentMasks.LongRun;
[long_run_segments, num_of_long_runs] = SegmentVector(long_run);

drop_center = movie_features.drop{1};

fp = track.filteredPath(:,2:3);
weathervane = [];
for segment_ix = 1:num_of_long_runs
    if ~long_run_segments(segment_ix,3)
        continue
    end
    
    start_ix = long_run_segments(segment_ix,1);
    stop_ix = long_run_segments(segment_ix,2);
    
    for ix = start_ix:stop_ix-1
        curr_point = fp(ix,:);
        % Look ahead to 1 mm
        for next_ix = ix+1:stop_ix
            next_point = fp(next_ix,:);
            if Dist(curr_point, next_point) > movie_features.pixels_per_mm
                weathervane(end+1,:) = [ ...
                    track.bearings(ix) ...
                    GetAngle(track.velocityVectors(ix,:), next_point - curr_point) ...
                    track.filteredPath(ix,1) ...
                    Dist(curr_point, drop_center) ...
                ];
                break
            end
        end
    end
end

end
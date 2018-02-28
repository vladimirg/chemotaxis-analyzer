function weathervane = Weathervane(track, movie_features)

segmented_track = SegmentTrackPath(track, movie_features.fps);
long_run = bitand(segmented_track, SegmentMasks.LongRun) == SegmentMasks.LongRun;
[long_run_segments, num_of_long_runs] = SegmentVector(long_run);

drop_center = movie_features.drop{1};

% We trust the bearings since we smoothed the path, but we still want to
% look ahead 1 mm for the psi, as Iino did.

weathervane = [];
for segment_ix = 1:num_of_long_runs
    if ~long_run_segments(segment_ix,3)
        continue
    end
    
    start_ix = long_run_segments(segment_ix,1);
    stop_ix = long_run_segments(segment_ix,2);
    
    for ix = start_ix:stop_ix-1
        curr_point = track.filteredPath(ix,2:3);
        % Look ahead to 1 mm
        for next_ix = ix+1:stop_ix
            next_point = track.filteredPath(next_ix,2:3);
            if Dist(curr_point, next_point) > movie_features.pixels_per_mm
                weathervane(end+1,:) = [ ...
                    track.bearings(ix) ...
                    GetAngle(track.velocityVecotrs(ix,:), next_point - curr_point) ...
                    track.filteredPath(ix,1) ...
                    Dist(curr_point, drop_center) ...
                ];
                break
            end
        end
    end
end

% figure;
% bin_edges = -180:20:180;
% bearing_bins = discretize(weathervane(:,1), bin_edges);
% avg_per_bin = zeros(length(bin_edges)-1, 1);
% for bin_ix = 1:length(avg_per_bin)
%     avg_per_bin(bin_ix) = mean(weathervane(bearing_bins == bin_ix,2));
% end
% bar(bin_edges(1:end-1), avg_per_bin);

end
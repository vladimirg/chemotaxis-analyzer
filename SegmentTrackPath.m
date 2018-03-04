function segmentedPath = SegmentTrackPath(track, fps)
% Output is relative to velocity vector origins.
% If there are N steps, then the output length will be N-1.
% Segments each point into a path type (possibly more than one), as
% specified in SegmentMasks.
% The output is a bit vector, as it is possible for a point to be
% classified into a number of path types (e.g., a point can be part of
% short run within a pirouette, but not be a sharp turn).

% TODO: verify this is consistent in each and every strain, and between
% different experiments.
short_run_duration_in_secs = 6.05; % According to Pierece-Shimamura.
short_run_duration_in_frames = short_run_duration_in_secs * fps;

if isfield(track, 'filteredPath')
    fp = track.filteredPath;
else
    fp = track.path;
end
[steps, ~] = size(fp);
segmentedPath = FindSharpTurns(track, fps);

[segments, segment_count] = SegmentVector(bitand(segmentedPath, SegmentMasks.SharpTurn) > 0);

for segment_ix = 1:segment_count
    segment_start = segments(segment_ix,1);
    segment_stop = segments(segment_ix,2);
    is_turn = segments(segment_ix,3);
    
    if is_turn
        continue
    end
    
    % It's a run - is it short or long?
    first_frame = fp(segment_start, 1);
    last_frame = fp(segment_stop, 1);
    % TODO: < or <= ?
    if last_frame - first_frame <= short_run_duration_in_frames
        bitmask = SegmentMasks.ShortRun;
    else
        bitmask = SegmentMasks.LongRun;
    end
    
    segmentedPath(segment_start:segment_stop) = bitor(...
        segmentedPath(segment_start:segment_stop),...
        bitmask);
end

[long_run_segments, long_run_count] = SegmentVector(bitand(segmentedPath, SegmentMasks.LongRun) > 0);

% A pirouette is any segment that contains at least one short run and no
% long runs.
% two_secs = ceil(fps * 2);
% points_count = 6;
for segment_ix = 1:long_run_count
    segment_start = long_run_segments(segment_ix,1);
    segment_stop = long_run_segments(segment_ix,2);
    is_long_run = long_run_segments(segment_ix,3);
    
    if ~is_long_run
        segmentedPath(segment_start:segment_stop) = bitor(...
            segmentedPath(segment_start:segment_stop),...
            SegmentMasks.Pirouette);
%         % The pirouette is long enough, now check the angles before and
%         % after.
%         can_compute_angle = false;
%         if segment_start - points_count >= 1 && segment_stop + points_count + 1 <= steps
%             before_origin = fp(segment_start - points_count,2:3);
%             after_origin = fp(segment_stop+1,2:3);
%             angle_diff = GetAngle(fp(segment_start,2:3) - before_origin, fp(segment_stop+1+points_count,2:3) - after_origin);
%             can_compute_angle = true;
%         end
%         
%         pirouette_length = (segment_stop - segment_start + 1);
%         
% %         if abs(angle_diff) <= 30
% %             continue
% %         end
% 
%         if pirouette_length <= two_secs && can_compute_angle && abs(angle_diff) <= 30
%             bitmask = SegmentMasks.Reversal;
% %         elseif pirouette_length <= 2
% %             % If we have 2 frames or less, it's probably an omega turn.
% %             % TODO: not really how this works.
% %             bitmask = SegmentMasks.OmegaTurn;
%         else
%             bitmask = SegmentMasks.Pirouette;
%         end
%         
%         segmentedPath(segment_start:segment_stop) = bitor(...
%             segmentedPath(segment_start:segment_stop),...
%             bitmask);
    end
end

end
function tracking_info = AnalyzeTrackingInfo(filename, movie_features)

load(filename, 'tracks', 'tracker');

plate_center = movie_features.plate{1};
plate_radius = movie_features.plate{2};
drop_center = movie_features.drop{1};
drop_radius = movie_features.drop{2};

filtered_tracks = [];

for track_ix = 1:length(tracks)
    % Filtering criteria:
    % 1) Tracks must be entirely within the plate, but outside the drop.
    % 2) At least some of the track should fall within the the frame range
    % specified in the frames file -- NOT ENFORCED HERE ANYMORE!!!
    % 3) The median displacement of the track from the origin must be
    % larger than 1 mm.
    
    track = tracks(track_ix);
    path = track.path;
    
%     if path(1,1) > stop_frame || path(end,1) < start_frame
%         continue
%     end
    
    if any(GetDistanceFromPoint(path(:,2:3), plate_center) > plate_radius) || ...
        any(GetDistanceFromPoint(path(:,2:3), drop_center) < drop_radius)
        continue;
    end
    
    origin = path(1,2:3);
    displacement_from_origin = cellfun(@(point) Dist(origin, point(2:3)), num2cell(path, 2));
    median_displacement = median(displacement_from_origin);
    if median_displacement <= movie_features.pixels_per_mm
        continue
    end
    
    path = [path(:,1) smooth(path(:,2)) smooth(path(:,3)) path(:,4)];
    
    % TODO: don't filter for now, see how that works out.
%     filtered_path = FilterPath(track);
%     filtered_step_sizes = StepSizes(filtered_path);
    track.filteredPath = path;
    track.filteredStepSizes = StepSizes(path);
    velocity_vectors = GetVelocityVectors(path);
    track.velocityVectors = velocity_vectors;
    track.bearings = GetBearingsToTargetPoint(path, drop_center);
    track.absoluteBearings = GetBearingsToVector(path, [1 0]);
    track.anglesBetweenSteps = ...
        GetAnglesBetweenSteps(velocity_vectors);
    
    % Sometimes worms pass through a stain on the plate, even after
    % filtering. The tracker thinks the worm remained in place, and this
    % can create hundreds and thousands of near-stationary points. We
    % filter them out by using a moving window - the values were chosen
    % empirically.
    jump_ix = 50;
    fp = path(:,2:3);
    [num_of_points, ~] = size(fp);
    bad_region = false;
    for p_ix = jump_ix+1:jump_ix:(num_of_points-jump_ix)
        origin = fp(p_ix,:);
        start_ix = p_ix - jump_ix;
        stop_ix = p_ix + jump_ix;
        dists_from_origin = arrayfun(@(a,b) Dist(origin, [a b]), fp(start_ix:stop_ix,1), fp(start_ix:stop_ix,2));
        if median(dists_from_origin) <= 0.5
            bad_region = true;
            break;
        end  
    end
    if bad_region
        continue
    end
    
    % Worms can't move more than 1 mm per second (this estimate is very
    % conservative):
    if any(track.filteredStepSizes > (movie_features.pixels_per_mm / movie_features.fps))
        continue
    end
    
    if isempty(filtered_tracks)
        filtered_tracks = track;
    else
        filtered_tracks(end+1) = track;
    end
end

% The video file itself will have an .avi extension, which we have to
% replace with .mat.
[file_folder, video_filename, ~] = fileparts(tracker.videoFileName);
video_comments = fullfile(file_folder, [video_filename(1:end-4) '.mat']);
load(video_comments);

% Add filename, concentration, date, strain.
% tracking_info = struct();
tracking_info = struct(...
    'tracks', filtered_tracks, ...
    'tracker', tracker, ...
    'filename', tracker.name, ...
    'comments', parameters ...
);

end

function filtered_path = FilterPath(track, fps, px_per_mm)

max_velocity_in_mm = 1;
max_velocity_in_px = max_velocity_in_mm * px_per_mm / fps;

min_dist_in_px = 0.01; % TODO: how do we handle standing worms?
stepSizes = StepSizes(track.path);
% TODO: there is a problem here of jumps - non-transient jumps will not
% be discarded!
filtered_path = track.path([true; (stepSizes <= max_velocity_in_px) & (stepSizes > min_dist_in_px)],:);

end

function stepSizes = StepSizes(path)

stepSizes = sqrt((path(1:end-1,2) - path(2:end,2)) .^ 2 +...
    (path(1:end-1,3) - path(2:end,3)) .^ 2);

end

function speed_vecs = GetVelocityVectors(path)

speed_vecs = path(2:end,2:3) - path(1:end-1,2:3);

end

function bearings = GetBearingsToTargetPoint(path, target_point)

origins = path(1:end-1,2:3);
[num_of_origins, ~] = size(origins);
spatial_vectors = path(2:end,2:3) - origins;
vectors_to_target = repmat(target_point, length(origins), 1) - origins;

bearings = zeros(length(origins), 1);
for ix = 1:num_of_origins
    bearings(ix) = GetAngle(vectors_to_target(ix,:), spatial_vectors(ix,:));
end

end

function bearings = GetBearingsToVector(path, ref_vector)

origins = path(1:end-1,2:3);
velocity_vectors = path(2:end,2:3) - origins;
[num_of_velocity_vector, ~] = size(origins);

bearings = zeros(num_of_velocity_vector, 1);
for ix = 1:num_of_velocity_vector
    bearings(ix) = GetAngle(ref_vector, velocity_vectors(ix,:));
end

end

function dists = GetDistanceFromPoint(path, point)

[path_length, ~] = size(path);
dists = zeros(path_length, 1);
for p_ix = 1:path_length
    dists(p_ix) = Dist(path(p_ix,:), point);
end

end
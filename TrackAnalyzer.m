function TrackAnalyzer(in_folder, movie_features_file, out_folder)

min_median_displacement = 10;
min_track_duration = 60;

mat_files = GetTrackingInfoFiles(in_folder);
load(movie_features_file);
all_features = features; % For parfor

frames_file = 'F:\tmp\test tracks\frames.xlsx';
[~, ~, frames_raw] = xlsread(frames_file, 'Comments');

% analyzed_data = struct();

header = {'Filename', 'Strain', 'Odorant concentration', 'Filtered tracks',...
    'Mean length', 'Median length',...
    'Mean steps', 'Median steps', ...
    'Plate radius', 'Tracks >= 2/3 R', 'Tracks >= R'};
output = cell(length(mat_files) + 1, length(header));
output(1,:) = header;
for ix = 1:length(mat_files)
    tracking_info_file = mat_files{ix};
    
%     if ~(contains(tracking_info_file, 'Dec'))
%         continue
%     end
    
    frames_row = frames_raw(cellfun(@(n) contains(tracking_info_file, num2str(n)), frames_raw(:,1)) & ...
        (cellfun(@ischar, frames_raw(:,2)) |...
         ~cellfun(@(n) (~ischar(n) && isnan(n)), frames_raw(:,2))),...
        2:3);
    start_frame = frames_row{1};
    stop_frame = frames_row{2};
    if ischar(stop_frame)
        frame_split = strsplit(stop_frame, '-');
        stop_frame = str2num(frame_split{1});
    end
    
    movie_features = all_features(cellfun(@(name) contains(tracking_info_file, name(1:end-5)), {all_features.name}));
    tracking_info = AnalyzeTrackingInfo(tracking_info_file, movie_features, start_frame, stop_frame);
    
%     plate_center = movie_features.plate{1};
    plate_radius = movie_features.plate{2};
    
    tracks = tracking_info.tracks;
%     path_centroids = cellfun(@(p) mean([[0 0]; p(:,2:3)], 1), {tracks.filteredPath}, 'UniformOutput', false);
%     
%     filtered_tracks = tracks([tracks.medianDisplacement] >= min_median_displacement & ...
%         cellfun(@length, {tracks.filteredStepSizes}) >= min_track_duration & ...
%         cellfun(@(c) sqrt((c(1)-plate_center(1))^2 + (c(2)-plate_center(2))^2) <= plate_radius, path_centroids));
    track_lengths = cellfun(@sum, {tracks.filteredStepSizes});
    track_durations = cellfun(@length, {tracks.filteredStepSizes});
%     
%     drop_center = movie_features.drop{1};
%     for track_ix = 1:length(filtered_tracks)
%         track = filtered_tracks(track_ix);
%         velocity_vectors = GetVelocityVectors(track.filteredPath);
%         filtered_tracks(track_ix).velocityVecotrs = velocity_vectors;
%         filtered_tracks(track_ix).bearings = GetBearings(track.filteredPath, drop_center);
% %         filtered_tracks(track_ix).anglesFromDrop = ...
% %             GetAngles(track.filteredPath, drop_center);
%         filtered_tracks(track_ix).anglesBetweenSteps = ...
%             GetAnglesBetweenSteps(velocity_vectors);
%     end
    
%     [sorted_track_lengths, sorted_length_ixes] = sort(track_lengths);
%     sorted_tracks = filtered_tracks(sorted_length_ixes);
    
%     analyzed_data(ix).name = tracking_info.tracker.name;
%     analyzed_data(ix).tracks = filtered_tracks;
%     analyzed_data(ix).movieFeautres = movie_features
    name = tracking_info.tracker.name;
    comments = tracking_info.comments;
    tracker = tracking_info.tracker;
    save(fullfile(out_folder, [name '.analyzed_tracks.mat']), 'name', 'tracks', 'movie_features', 'comments', 'tracker');
    
    output(ix+1,:) = {tracking_info.tracker.name, tracking_info.comments.Strain, ...
        tracking_info.comments.OdorantConcentration, length(tracks),...
        mean(track_lengths), median(track_lengths),...
        mean(track_durations), median(track_durations), ...
        plate_radius, sum(track_lengths >= (plate_radius * 2 / 3)), ...
        sum(track_lengths >= plate_radius)};
    
    fprintf('Finished %d / %d\n', ix, length(mat_files));
end

xlswrite(fullfile(out_folder, 'track_analyzer_output.xlsx'), output);

end

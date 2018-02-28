function FeatureExtractor(in_folder, out_folder)

mat_files = GetTrackingInfoFiles(in_folder);

feature_file = fullfile(out_folder, 'MovieFeatures.mat');
if exist(feature_file, 'file')
    load(feature_file);
else
    features = struct();
end

for ix = 1:length(mat_files)
    tracking_info_file = mat_files{ix};
    load(tracking_info_file);
    
    % The video file itself will have an .avi extension, which we have to
    % replace with .mat.
    [file_folder, video_filename, ~] = fileparts(tracker.videoFileName);
    video_comments = fullfile(file_folder, [video_filename(1:end-4) '.mat']);
    
    existing_feature = cellfun(@(s) strcmp(s, tracker.name), {features.name});
    if any(existing_feature)
        % If FPS doesn't exist:
        feature_ix = find(existing_feature);
        if ~isfield(features, 'fps') || isempty(features(feature_ix).fps)
            load(video_comments);
            fps = tracker.numberOfFrames / parameters.Duration;
            
            features(feature_ix).fps = fps;
            features(feature_ix).pixels_per_mm = ...
                features(feature_ix).plate{2} * 2 / 90;
        end
        
        continue
    end
    
    % Load first image and show it, for the feature extraction to begin.
    video_filename = tracker.videoFileName;
    video = VideoReader(video_filename);
    first_frame = readFrame(video);
    
    if contains(tracking_info_file, '04-Feb-2018-18.42.52-Mic1-N2_-2_A')
        % I played with the zoom at the beginning of the video, and it
        % stabilized at the 300-th frame.
        video = VideoReader(video_filename);
        first_frame = read(video, 300);
    end
    
    imshow(first_frame);
    
    features(end+1).name = tracker.name;
    
    % We want two things - the plate and its center.
    [plate_x, plate_y] = getpts;
    [drop_x, drop_y] = getpts;
    close;
    
    [R, C, ~] = ExactMinBoundCircle(horzcat(plate_x, plate_y));
    features(end).plate = {C R};
    
    [R, C, ~] = ExactMinBoundCircle(horzcat(drop_x, drop_y));
    features(end).drop = {C R};
    
    load(video_comments);
    fps = tracker.numberOfFrames / parameters.Duration;
    features(end).fps = fps;
    features(end).pixels_per_mm = ...
        features(end).plate{2} * 2 / 90;
end

save(feature_file, 'features');

end
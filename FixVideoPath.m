function FixVideoPath(tracking_info_folder, video_folder)

for track_info_file_cell = GetTrackingInfoFiles(tracking_info_folder)
    track_info_file = track_info_file_cell{1};
    [base_folder, ~, ~] = fileparts(track_info_file);
    backup_file = fullfile(base_folder, sprintf('TrackingInfo.%s.bak.mat', datestr(datetime, 'yyyy-mm-dd_HH-MM')));
    
    load(track_info_file, 'tracker', 'tracks');
    [~, name, ext] = fileparts(tracker.videoFileName);
    new_video_filename = fullfile(video_folder, [name ext]);
    if ~exist(new_video_filename, 'file')
        fprintf('%s new video file %s does not exist - ignoring.\n', track_info_file, new_video_filename);
        continue;
    end
    
    copyfile(track_info_file, backup_file);
    tracker.videoFileName = new_video_filename;
    tracker.initialize;
    save(track_info_file, 'tracker', 'tracks');
end

end
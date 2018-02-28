function GoOverFilteredTracks
% TODO: actually we go over analyzed tracks, we filter them ourselves.

% mat_files = GetTrackingInfoFiles('F:\tmp\test tracks smoothed');
% out_folder = 'F:\tmp\test tracks smoothed\figs';
mat_files = GetTrackingInfoFiles('F:\tmp\test tracks');
out_folder = 'F:\tmp\test tracks\figs-feb';
frames_file = 'F:\tmp\test tracks\frames.xlsx';
[~, ~, frames_raw] = xlsread(frames_file, 'Comments');

header = {'Filename', 'Strain', 'Odorant concentration', 'Date',...
    'Total tracks', 'Filtered tracks', 'Mean length', 'Median length',...
    'Total steps', 'Long runs', 'Pirrouetes', 'Max. worms', '% in long runs', ...
    '% in pirouttes'};%, 'Chemotaxis coeff.'};
output = cell(length(mat_files) + 1, length(header));
output(1,:) = header;

header_chemo = {'Filename', 'Strain', 'Odorant concentration', 'Date', 'Max. worms', ...
    'First arrival', 'First 10%', 'First 20%', 'First 40%', 'First 80%', ...
    'Arrival at 5', 'Arrival at 10', 'Arrival at 15', 'Arrival at 20'...
    };
output_chemo = cell(length(mat_files) + 1, length(header_chemo));
output_chemo(1,:) = header_chemo;

combined_data = containers.Map;

% Which figures do we want to force?
force_bearings_fig = true;
force_polygon_dynamics_fig = false;
force_att_vec_field = false;
force_distances_fig = false;
force_turn_improvement = false;
force_top_10 = false;
force_top_10_segmented = false;
force_weathervane_fig = false;

for ix = 1:length(mat_files)
    output_steps = 0;
    output_pirouettes = 0;
    output_pirouette_steps = 0;
    output_kinks = 0;
    output_long_runs = 0;
    output_long_run_steps = 0;
    
    tracking_info_file = mat_files{ix};
    if ~contains(tracking_info_file, 'Feb')%'16.53.12');%'17.28.42')
        continue
    end
    load(tracking_info_file)

    frames_row = frames_raw(cellfun(@(n) contains(tracker.name, num2str(n)), frames_raw(:,1)) & ...
        (cellfun(@ischar, frames_raw(:,2)) |...
         ~cellfun(@(n) (~ischar(n) && isnan(n)), frames_raw(:,2))),...
        2:3);
    start_frame = frames_row{1};
    stop_frame = frames_row{2};
    if ischar(stop_frame)
        frame_split = strsplit(stop_frame, '-');
        stop_frame = str2num(frame_split{1});
    end
    total_frames = stop_frame - start_frame + 1;
        
    fps = movie_features.fps;
    
    % TODO: what would be a good number here?
    min_track_duration = ceil(30 * fps);
    
    % Filter the raw tracks into something sensible.
    plate_center = movie_features.plate{1};
    plate_radius = movie_features.plate{2};
    drop_center = movie_features.drop{1};
    drop_radius = movie_features.drop{2};
%     path_centroids = cellfun(@(p) mean([[0 0]; p(:,2:3)], 1), {tracks.filteredPath}, 'UniformOutput', false);
    
    px_per_mm = movie_features.pixels_per_mm;
    
    filtered_tracks = tracks(...
        cellfun(@(p) p(end,1) - p(1,1), {tracks.filteredPath}) >= min_track_duration & ...
        cellfun(@(p) Dist(drop_center, p(1,2:3)), {tracks.filteredPath}) > 1/4*plate_radius & ...
        cellfun(@(p) p(end,1) > start_frame, {tracks.filteredPath}) & ...
        cellfun(@(p) p(1,1) < stop_frame, {tracks.filteredPath})...
        );

    % Strain - concentration - date
    strain = comments.Strain;
    odorant_concentration = comments.OdorantConcentration(end-1:end);
    movie_date = comments.Date;
    label = sprintf('%s - %s - %s', strain, odorant_concentration, movie_date);
    filename_label = sprintf('%s - %s - %s', strain, odorant_concentration, name);
    
    track_lengths = cellfun(@sum, {filtered_tracks.filteredStepSizes});
    track_durations = cellfun(@length, {filtered_tracks.filteredStepSizes});
    
    [sorted_track_lengths, sorted_length_ixes] = sort(track_lengths);
    sorted_tracks = filtered_tracks(sorted_length_ixes);
    
%     longest_tracks = filtered_tracks(track_durations > 5*60*fps & ...
%         cellfun(@(p) Dist(drop_center, p(1,2:3)), {filtered_tracks.filteredPath}) > plate_radius * 4/5);
%     PlotTracks(longest_tracks, movie_features);
%     title(label);
%     savefig(fullfile(out_folder, [filename_label '.5_min_tracks.fig']));
%     close;

    bearings = [];
    bearing_distances = [];
    long_run_bearings = [];
    long_run_speeds = [];
    
    % Prepirouettes are a 3D matrix:
    % x - prepirouette #
    % y - step ix
    % z - bearing (1) or speed (2)
    prepirouette_length = round(fps * 60); % seconds
    prepirouettes = zeros(0, prepirouette_length, 2);
    
    % Look only at tracks the start location of which is beyond half a
    % radius:
    start_distances = cellfun(@(p) sqrt((plate_center(1) - p(1,2))^2 + (plate_center(2) - p(1,3))^2), {sorted_tracks.filteredPath});
%     sorted_tracks = sorted_tracks(start_distances < plate_radius*4/5 & start_distances > plate_radius*1/5 );
    sorted_tracks = sorted_tracks(start_distances > plate_radius * 1/5);
    
    for track_ix = 1:length(sorted_tracks)
        track = sorted_tracks(track_ix);
        fp = track.filteredPath(:,2:3);
        [num_of_steps, ~] = size(fp);
        segmented_path = SegmentTrackPath(track, fps);
        
        long_runs_all = bitand(segmented_path, SegmentMasks.LongRun) > 0 & ...
            bitand(segmented_path, SegmentMasks.Kink) == 0;
        long_run_bearings = [long_run_bearings track.bearings(long_runs_all)'];
        long_run_speeds = [long_run_speeds track.filteredStepSizes(long_runs_all)'];
        
        [pirouettes, num_of_pirouettes] = SegmentVector(bitand(segmented_path, SegmentMasks.Pirouette) > 0);
        % We only want to consider pirouettes that had a long run before
        % and after.
%         if num_of_pirouettes < 4 || (num_of_pirouettes == 3 && pirouttes(1,3) == 1)
%             continue
%         end
        
        output_steps = output_steps + length(track.filteredStepSizes);
        output_pirouettes = output_pirouettes + sum(pirouettes(:,3) > 0);
        output_pirouette_steps = output_pirouette_steps + sum(bitand(segmented_path, SegmentMasks.Pirouette) > 0);
        output_kinks = output_kinks + sum(bitand(segmented_path, SegmentMasks.Kink) > 0);
        output_long_runs = output_long_runs + sum(pirouettes(:,3) == 0);
        output_long_run_steps = output_long_run_steps + sum(cellfun(@(p) (p(3) == 0) * (p(2) - p(1) + 1), mat2cell(pirouettes, ones(num_of_pirouettes, 1))));
        
        for pix = 2:num_of_pirouettes
            if pirouettes(pix,3) == 0
                continue
            end
            
            pirouette_start = pirouettes(pix,1);
            pirouette_end = pirouettes(pix,2);
            
            long_run_start = pirouettes(pix-1,1);
            long_run_end = pirouettes(pix-1,2);
            if long_run_end - long_run_start >= prepirouette_length
                prepirouette = ones(prepirouette_length, 2) * nan;
                prepirouette(:,1) = track.bearings(long_run_end-prepirouette_length+1:long_run_end);
                prepirouette(:,2) = track.filteredStepSizes(long_run_end-prepirouette_length+1:long_run_end) / px_per_mm * fps;
                prepirouettes(end+1,:,:) = prepirouette;
            end
            
            % New method
            % TODO: document
            points_count = 6;
            before_origin = fp(pirouette_start - points_count,:);
            bearing_before_pir = GetAngle(drop_center - before_origin, fp(pirouette_start,:) - before_origin);
            
            if pix < num_of_pirouettes
                after_origin = fp(pirouette_end+1,:);
                bearing_after_pir = GetAngle(drop_center - after_origin, fp(pirouette_end+1+points_count,:) - after_origin);
                bearing_distances = [bearing_distances; Dist(drop_center, fp(pirouette_start,:)) Dist(drop_center, after_origin)];
            else
                bearing_after_pir = nan;
            end
            
            bearings = [bearings; [bearing_before_pir bearing_after_pir]];
        end
    end
    
    matched_bearings = bearings(~isnan(bearings(:,2)),:);
    
    bin_edges = -180:20:180;
    
    % What are the speed averages, relative to bearing?
    long_run_speeds = long_run_speeds / px_per_mm * fps;
    binned_bearings = discretize(long_run_bearings, bin_edges);
    binned_speeds = [];
    for bin_ix = 1:length(bin_edges) - 1
        binned_speeds(end+1) = mean(long_run_speeds(binned_bearings == bin_ix));
    end
    
    abs_speed_fig = [fullfile(out_folder, filename_label) '.speeds_vs_bearing.fig'];
    if ~exist(abs_speed_fig, 'file')
        figure;
        bar(bin_edges(1:end-1), binned_speeds);
        title(label);
        savefig(abs_speed_fig);
        close;
    end
    
    normed_bearings = GetConditionalTurnProbabilities(long_run_bearings, bearings(:,1), bin_edges);
    
    distances_per_frame = GetDistancesPerFrame(tracks, drop_center);
    tracked_worms_per_frame = sum(~isnan(distances_per_frame), 1);
    
    fig_filename = [fullfile(out_folder, filename_label) '.fig'];
    if force_bearings_fig || ~exist(fig_filename, 'file')
        fig_rows = 4;
        fig_cols = 2;
        
        hists = figure;
        
        subplot(fig_rows, fig_cols, 1);
        histogram(long_run_bearings, bin_edges);
        title('# steps vs. B');

        subplot(fig_rows, fig_cols, 2);
        histogram(bearings(:,1), bin_edges);
        title('# pir vs. B_{before}');

        subplot(fig_rows, fig_cols, 3);
        bar(bin_edges(1:end-1), normed_bearings);
        title('P(pir|B_{before})');

        subplot(fig_rows, fig_cols, 4);
        histogram(matched_bearings(:,2), bin_edges);
        title('# pir vs. B_{after}');

        binned_bearings = discretize(matched_bearings(:,1), bin_edges);
        deltas = arrayfun(@(b, a) GetAngle([cosd(b) sind(b)], [cosd(a) sind(a)]), ...
            matched_bearings(:,1), matched_bearings(:,2));
        abs_delta_bearings = [];
        for bin_ix = 1:length(bin_edges) - 1
            abs_delta_bearings(end+1) = mean(abs(deltas(binned_bearings == bin_ix)));
        end
        
        subplot(fig_rows, fig_cols, 5);
        histogram(deltas, bin_edges);
        title('# pir vs. dB');
        
        subplot(fig_rows, fig_cols, 6);
        bar(bin_edges(1:end-1), abs_delta_bearings);
        title('Abs mean dB vs. B_{before}');
        
        ts = linspace(-prepirouette_length / fps, 0, prepirouette_length);
        
        subplot(fig_rows, fig_cols, 7);
        plot(ts, mean(abs(prepirouettes(:,:,1))));
        hold on;
        plot(ts, std(abs(prepirouettes(:,:,1))));
%         lgd = legend('B', 'std dev');
%         lgd.Box = 'off';
%         lgd.Location = 'northwest';
        xlim([-60 0]);
        title('Prepirouette B vs. t (sec)');
        
        subplot(fig_rows, fig_cols, 8);
        plot(ts, mean(prepirouettes(:,:,2)));
        xlim([-60 0]);
        title('Prepirouette speed vs. t (sec)');

        mtit(hists, label, 'xoff', 0, 'yoff', 0.03, 'zoff', 0);
        savefig(hists, fig_filename);
        close;
    end

    chem_coeff = nan;
    % TODO: a mess!
    pol_dyn_fig_name = fullfile(out_folder, [filename_label '.norm_polygon_dynamics.fig']);
    if force_polygon_dynamics_fig || ~exist(pol_dyn_fig_name, 'file') 
        drop_center = movie_features.drop{1};
        drop_radius = movie_features.drop{2};
        worms_in_region = polygonDynamics(tracker, 1, 'Polygon', circleToPolygon([drop_center drop_radius*2.5], 8),...
            'Tracks', tracks);
        title([label ' - ' name]);
%         savefig(pol_dyn_fig_name);
        
        % Find first and max indexes:
%         max_tracked_worms = max(sorted_tracks_tracked_worms);
%         max_worms = max(worms_in_region);
%         worms_reached_ixes = find(worms_in_region >= max_worms*0.1);
%         first_ix = worms_reached_ixes(1);
%         worms_max_ixes = find(worms_in_region >= max_worms*0.9);
%         last_ix = worms_max_ixes(1);
%         worms_slice = worms_in_region(first_ix:last_ix) / max_tracked_worms;
        
        close;
        
%         chemo_velocity_fig_name = fullfile(out_folder, [filename_label '.chemo_velocity.fig']);
%         figure;
%         
%         x_axis_values = 1:(last_ix-first_ix+1);
%         plot(x_axis_values, worms_slice, 'o');
%         [coeffs, ~] = polyfit(x_axis_values, worms_slice, 1);
%         chem_coeff = coeffs(1);
%         y_axis_intersection = coeffs(2);
%         
%         hold on;
%         plot(x_axis_values, (x_axis_values * chem_coeff) + y_axis_intersection);
%         
%         title(sprintf('%s - %f', filename_label, chem_coeff));
%         
%         savefig(chemo_velocity_fig_name);
%         close;
%         
%         figure;
        
        normed_worms_in_region = worms_in_region(start_frame:end) / max(tracked_worms_per_frame);
        ts = linspace(0, length(normed_worms_in_region)/fps, length(normed_worms_in_region)) / 60;
        stop_time = (stop_frame-start_frame) / fps / 60;
        
        first_worm_time = find(normed_worms_in_region > 0) / fps / 60;
        first_worm_time = first_worm_time(1);

        figure;
        hold on;
        plot(ts, normed_worms_in_region);
        plot([stop_time stop_time], [0, 1.1], 'k--')
        title(label);
        
        new_ticks = unique(sort([first_worm_time, stop_time, get(gca, 'XTick')]));
        new_labels = arrayfun(@num2str, new_ticks, 'UniformOutput', false);
        new_labels{new_ticks == first_worm_time} = sprintf('%.1f', first_worm_time);
        new_labels{new_ticks == stop_time} = sprintf('%.1f', stop_time);
        set(gca, 'XTick', new_ticks);
        set(gca, 'XTickLabel', new_labels);
        ylim([0 1.1]);
        
        savefig(fullfile(out_folder, [filename_label '.norm_polygon_dynamics.fig']));
        close;
        
%         header_chemo = {'Filename', 'Strain', 'Odorant concentration', 'Max. worms', ...
%             'First arrival', 'First 10%', 'First 20%', 'First 40%', 'First 80%', ...
%             'Arrival at 5', 'Arrival at 10', 'Arrival at 15', 'Arrival at 20'...
%             };
%         output_chemo = cell(length(mat_files) + 1, length(header_chemo));
%         output_chemo(1,:) = header_chemo;

%         output(ix+1,:) = {tracker.name, comments.Strain, ...
%             comments.OdorantConcentration, movie_date,...
%             length(tracks), length(sorted_tracks), mean(sorted_track_lengths), median(sorted_track_lengths),...
%             output_steps, output_long_runs, output_pirouettes,...
%             max(tracked_worms_per_frame), output_long_run_steps / output_steps,...
%             output_pirouette_steps / output_steps, chem_coeff...
%         };
        
        first_arrivals = {};
        for arrival_threshold = [0.1 0.2 0.4 0.8]
            arrival_time = nan;
            arrival_ix = min(find(normed_worms_in_region >= arrival_threshold));
            if ~isempty(arrival_ix)
                arrival_time = datestr(duration(0, ts(arrival_ix), 0), 'MM:SS');
            end
            first_arrivals{end+1} = arrival_time;
        end
        
        arrived_percent = [];
        for arrival_time = [5 10 15 20]
            arrived_percent(end+1) = normed_worms_in_region(round(arrival_time * 60 * fps));
        end
        
        output_chemo(ix+1,:) = {tracker.name, comments.Strain, comments.OdorantConcentration, ...
            movie_date, max(tracked_worms_per_frame), datestr(duration(0, first_worm_time, 0), 'MM:SS'), ...
            first_arrivals{1}, first_arrivals{2}, first_arrivals{3}, first_arrivals{4}, ...
            arrived_percent(1), arrived_percent(2), arrived_percent(3), arrived_percent(4) ...
        };
    end
    
    att_vec_field_fig_name = fullfile(out_folder, [filename_label '.attraction_vector_field.fig']);
    if force_att_vec_field || ~exist(att_vec_field_fig_name, 'file') 
        attractionVectorField(tracker, 0, 'Tracks', sorted_tracks');
        title([label ' - ' name]);
        savefig(att_vec_field_fig_name);
        close;
    end
    
    long_tracks = sorted_tracks(cellfun(@(fp) Dist(drop_center, fp(1,2:3)) >= 4/5*plate_radius, {sorted_tracks.filteredPath}) & ...
        cellfun(@(fp) Dist(drop_center, fp(end,2:3)) <= 1/5*plate_radius, {sorted_tracks.filteredPath}));
    
    mean_velocities = [];
    for track_ix = 1:length(sorted_tracks)
        track = sorted_tracks(track_ix);
        if length(track.filteredStepSizes) < 60*fps || ...
            Dist(drop_center, track.filteredPath(1,2:3)) < plate_radius / 2
            continue;
        end
        
%         mean_velocities(end+1) = -(Dist(drop_center, track.filteredPath(end,2:3)) - ...
%             Dist(drop_center, track.filteredPath(1,2:3))) / px_per_mm / ...
%             ((track.filteredPath(end,1) - track.filteredPath(1,1) + 1) / fps);
        mean_velocities(end+1) = GetAvgVelocity(track, drop_center, px_per_mm, fps);
    end
    
%     figure;
%     histogram(mean_velocities, -0.1:0.005:0.1);
%     title(label);
%     savefig(fullfile(out_folder, [filename_label '.velocity_histogram.fig']));
%     close;
    
    med_dist_plot_name = fullfile(out_folder, [filename_label '.all_distances.fig']);
    if force_distances_fig || ~exist(med_dist_plot_name, 'file')
        [distances, start_frame, stop_frame] = GetDistancesPerFrame(tracks, drop_center);
        tracked_worms = sum(~isnan(distances), 1);
        
        % Report the distance data.
        figure;
        hold on;
        title(label);
        mean_dists = mean(distances, 1, 'omitnan') / plate_radius;
        median_dists = median(distances, 1, 'omitnan') / plate_radius;
        std_dev = std(distances, 0, 1, 'omitnan') / plate_radius;
        tracked_worms = tracked_worms / max(tracked_worms); % Normalize
        xs = start_frame:stop_frame;
        plot(xs, mean_dists, 'b');
        plot(xs, median_dists, 'r');
        plot(xs, std_dev, 'g');
        plot(xs, tracked_worms, 'm');

        new_ticks = unique(sort([start_frame, stop_frame, get(gca, 'XTick')]));
        new_labels = arrayfun(@num2str, new_ticks, 'UniformOutput', false);
        new_labels{new_ticks == start_frame} = sprintf('%d S', start_frame);
        new_labels{new_ticks == stop_frame} = sprintf('%d E', stop_frame);
        set(gca, 'XTick', new_ticks);
        set(gca, 'XTickLabel', new_labels);
        ylim([0 1]);
        
        plot(1:tracker.numberOfFrames, ones(tracker.numberOfFrames, 1) * drop_radius / plate_radius, '--k');
        legend('mean', 'median', 'std dev', '# tracks', 'drop radius');
        
        savefig(med_dist_plot_name);
        close;
    end
    
    output(ix+1,:) = {tracker.name, comments.Strain, ...
        comments.OdorantConcentration, movie_date,...
        length(tracks), length(sorted_tracks), mean(sorted_track_lengths), median(sorted_track_lengths),...
        output_steps, output_long_runs, output_pirouettes,...
        max(tracked_worms_per_frame), output_long_run_steps / output_steps,...
        output_pirouette_steps / output_steps, ...%chem_coeff...
        };
    
    % Collect data per condition - strain and concentration.
    condition_key = sprintf('%s - %s', comments.Strain, comments.OdorantConcentration);
    if ~isKey(combined_data, condition_key)
        combined_data(condition_key) = struct('all_bearings', [], 'turn_bearings', [], ...
            'mean_velocities', [], 'weathervane', []);
    end
    
    all_weathervanes = [];
    for track_ix = 1:length(sorted_tracks)
        weathervane = Weathervane(sorted_tracks(track_ix), movie_features);
        if isempty(weathervane)
            continue
        end
        
        weathervane = weathervane( ...
            weathervane(:,3) - start_frame < fps * 60 * 15 & ...
            weathervane(:,4) <= 2 * plate_radius / 3 & ...
            weathervane(:,4) >= 1 * plate_radius / 3 ...
        );
        all_weathervanes = [all_weathervanes; Weathervane(sorted_tracks(track_ix), movie_features)];
    end
    
    prev_data = combined_data(condition_key);
    combined_data(condition_key) = struct(...
        'all_bearings', [prev_data.all_bearings long_run_bearings], ...
        'turn_bearings', [prev_data.turn_bearings; bearings], ...
        'mean_velocities', [prev_data.mean_velocities mean_velocities], ...
        'weathervane', [prev_data.weathervane; all_weathervanes]);
    
    turn_improvement_filename = fullfile(out_folder, [filename_label '.turn_improvement.fig']);
    if force_turn_improvement || ~exist(turn_improvement_filename, 'file')
        filtered_bearings = matched_bearings(abs(matched_bearings(:,1)) >= 90 ,:);
        theta_before = ToPolar(filtered_bearings(:,1));
        theta_after = ToPolar(filtered_bearings(:,2));
        
        polarhistogram(theta_before, 18);
        hold on;
        polarhistogram(theta_after, 18);
        set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');
        title(label);
        legend('Before', 'After');
        savefig(turn_improvement_filename);
        close;
    end
    
    top_10_tracks_filename = fullfile(out_folder, [filename_label '.top_10_tracks.fig']);
    if force_top_10 || ~exist(top_10_tracks_filename, 'file')
        top_10_tracks_fig = PlotTracks(sorted_tracks(end-9:end), movie_features);
        title(top_10_tracks_fig.CurrentAxes, label);
        savefig(top_10_tracks_fig, top_10_tracks_filename);
        close(top_10_tracks_fig);
    end
    
    top_10_segmented_tracks_filename = fullfile(out_folder, [filename_label '.segmented_top_10_tracks.fig']);
    if force_top_10_segmented || ~exist(top_10_segmented_tracks_filename, 'file')
        top_10_tracks_fig = PlotTracks(sorted_tracks(end-9:end), movie_features, 'SegmentWithFps', fps);
        title(top_10_tracks_fig.CurrentAxes, label);
        savefig(top_10_tracks_fig, top_10_segmented_tracks_filename);
        close(top_10_tracks_fig);
    end
    
    all_weathervane_figure = fullfile(out_folder, [filename_label '.all_weathervanes.fig']);
    if force_weathervane_fig || ~exist(all_weathervane_figure, 'file')
        w_bin_edges = -180:20:180;
        w_bearing_bins = discretize(all_weathervanes(:,1), w_bin_edges);
        w_avg_per_bin = zeros(length(w_bin_edges)-1, 1);
        for bin_ix = 1:length(w_avg_per_bin)
            w_avg_per_bin(bin_ix) = mean(all_weathervanes(w_bearing_bins == bin_ix,2));
        end
        
        figure;
        bar(w_bin_edges(1:end-1), w_avg_per_bin)
        title(label);

        savefig(all_weathervane_figure);
        close;
    end
    
    abs_speeds_filename = fullfile(out_folder, [filename_label '.abs_speeds.fig']);
    if ~exist(abs_speeds_filename, 'file')
        all_speeds = [];
        long_run_speeds = [];
        for track_ix = 1:length(sorted_tracks)
            track = sorted_tracks(track_ix);
            segmented_path = SegmentTrackPath(track, fps);
            [long_runs, num_of_segments] = SegmentVector(bitand(segmented_path, SegmentMasks.LongRun) > 0);

            all_speeds = [all_speeds; track.filteredStepSizes];
            for seg_ix = 1:num_of_segments
                if ~long_runs(seg_ix,3)
                    pir_speeds = track.filteredStepSizes(long_runs(seg_ix,1):long_runs(seg_ix,2))  / px_per_mm * fps;
                    if mean(pir_speeds) <= 0.01
                        1;
                    end
                    continue
                end

                long_run_speeds = [long_run_speeds; track.filteredStepSizes(long_runs(seg_ix,1):long_runs(seg_ix,2))];
            end
        end

        all_speeds = all_speeds / px_per_mm * fps;
        long_run_speeds = long_run_speeds / px_per_mm * fps;

        edges = 0:0.005:0.5;
        subplot(2, 1, 1);
        histogram(all_speeds, edges);
        title(label);

        subplot(2, 1, 2);
        histogram(long_run_speeds, edges);

        savefig(abs_speeds_filename);
        close;
    end
    
    fprintf('Finished %02d/%02d - %s\n', ix, length(mat_files), label);
%     break
end

return

xlswrite(fullfile(out_folder, 'movie_statistics.xlsx'), output);
xlswrite(fullfile(out_folder, 'chemotaxis_statistics.xlsx'), output_chemo);

% return
for condition_cell = keys(combined_data)
    condition = condition_cell{1};
    data = combined_data(condition);
    
    
    
    mean_velocities = data.mean_velocities;
    figure;
    histogram(mean_velocities, -0.1:0.005:0.1);
    title(condition);
    savefig(fullfile(out_folder, [condition '.combined.mean_velocities.fig']));
    close;
%     continue;



    all_weathervanes = data.weathervane;
    w_bin_edges = -180:20:180;
    w_bearing_bins = discretize(all_weathervanes(:,1), w_bin_edges);
    w_avg_per_bin = zeros(length(w_bin_edges)-1, 1);
    for bin_ix = 1:length(w_avg_per_bin)
        w_avg_per_bin(bin_ix) = mean(all_weathervanes(w_bearing_bins == bin_ix,2));
    end

    figure;
    bar(w_bin_edges(1:end-1), w_avg_per_bin)
    title(condition);

    savefig(fullfile(out_folder, [condition '.combined.weathervane.fig']));
    close;

    
    
    
    matched_bearings = data.turn_bearings(~isnan(data.turn_bearings(:,2)),:);
    filtered_bearings = matched_bearings(abs(matched_bearings(:,1)) >= 90 ,:);
    theta_before = ToPolar(filtered_bearings(:,1));
    theta_after = ToPolar(filtered_bearings(:,2));

    polarhistogram(theta_before, 18);
    hold on;
    polarhistogram(theta_after, 18);
    set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');
    title(condition);
    legend('Before', 'After');
    savefig(fullfile(out_folder, [condition '.combined.turn_improvement.fig']));
    close;
        
%     continue
    
    
    
    
    bin_edges = -180:20:180;
    
    conditional_probs = GetConditionalTurnProbabilities(data.all_bearings, data.turn_bearings(:,1), bin_edges);
    
    hists = figure;

    subplot(4, 1, 1);
    histogram(data.all_bearings, bin_edges);
    title('All bearings with turns');

    subplot(4, 1, 2);
    histogram(data.turn_bearings(:,1), bin_edges);
    title('Bearings before turn');

    subplot(4, 1, 3);
    bar(bin_edges(1:end-1), conditional_probs);
    title('Bearings before turns - conditional');

    subplot(4, 1, 4);
    histogram(data.turn_bearings(:,2), bin_edges);
    title('Bearings after');

    mtit(hists, condition, 'xoff', 0, 'yoff', 0.03, 'zoff', 0);
    savefig(hists, fullfile(out_folder, [condition '.combined.fig']));
    close;
end

end




function theta = ToPolar(angles)

angles = angles / 360 * 2 * pi;
% angles(angles < 0) = abs(angles(angles < 0)) + pi;

theta = angles;

end

function [mat_files] = GetTrackingInfoFiles(folder)

mat_files = {};

listing = dir(folder);

for ix = 1:numel(listing)
    entry = listing(ix);
    if entry.isdir && ~strcmp(entry.name, '.') && ~strcmp(entry.name, '..')
        mat_files = horzcat(mat_files, GetTrackingInfoFiles(fullfile(entry.folder, entry.name)));
    elseif contains(entry.name, '.analyzed_tracks.mat')
        mat_files{end+1} = fullfile(entry.folder, entry.name);
    end
end

end

function conditional_bins = GetConditionalTurnProbabilities(all_bearings, bearings_before_turn, bin_edges)

% Normalize the 'before' bearings to get probabilities of turns
before_bearings_bins = discretize(bearings_before_turn, bin_edges);
long_run_bearings_bins = discretize(all_bearings, bin_edges);

normed_bearings = zeros(length(bin_edges) - 1, 1);
for bin_ix = 1:length(normed_bearings)
    turns_at_bearing = sum(before_bearings_bins == bin_ix);
    steps_at_bearing = sum(long_run_bearings_bins == bin_ix);
    normed_bearings(bin_ix) = turns_at_bearing / steps_at_bearing;
end

conditional_bins = normed_bearings;

end

function [distances, first_frame, last_frame] = GetDistancesPerFrame(tracks, target_point)

first_frame = min(cellfun(@(fp) fp(1,1), {tracks.filteredPath}));
last_frame = max(cellfun(@(fp) fp(end,1), {tracks.filteredPath}));
num_of_frames = last_frame - first_frame + 1;
num_of_tracks = length(tracks);

distances = ones(num_of_tracks, num_of_frames) * nan;

for track_ix = 1:num_of_tracks
    track = tracks(track_ix);
    [num_of_points, ~] = size(track.filteredPath);
    for point_ix = 1:num_of_points
        point = track.filteredPath(point_ix,:);
        distances(track_ix, point(1) - first_frame + 1) = Dist(point(2:3), target_point);
    end
end

end

function avg_velocity = GetAvgVelocity(track, ref_point, px_per_mm, fps)

velocities = [];

segmented_track = SegmentTrackPath(track, fps);
[long_runs, num_of_long_runs] = SegmentVector(bitand(segmented_track, SegmentMasks.LongRun) > 0);

for segment_ix = 1:num_of_long_runs
    if ~long_runs(segment_ix,3)
        continue
    end
        
    for vec_ix = long_runs(segment_ix,1):long_runs(segment_ix,2)
    %     vel_vec = track.velocityVecotrs(vec_ix,:);
    %     target_vector = ref_point - track.filteredPath(vec_ix,2:3);
    %     proj_vec = dot(vel_vec,target_vector)/norm(target_vector)^2*target_vector;
    %     velocities(vec_ix) = sqrt( proj_vec(1)^2 + proj_vec(2)^2 );

    %     vel_vec_length = sqrt( vel_vec_length(1)^2 + vel_vec_length(2)^2 ) / px_per_mm;
        velocities(end+1) = (track.filteredStepSizes(vec_ix) / px_per_mm) * ...
            cos(deg2rad(track.bearings(vec_ix)));
    end
end

avg_velocity = mean(velocities);

end
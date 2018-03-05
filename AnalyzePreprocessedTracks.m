function AnalyzePreprocessedTracks(data_folder, out_folder, varargin)
% Analyze files with pre-processed tracks, and output results into a given
% folder.
% Optionally can accept a filter function for the the input files as a
% 'FilterByName' parameter, which accepts the name of the file, and returns
% true or false (an output of false means that the 

%% Parse optional arguments:
p = inputParser;
addParameter(p, 'FilterByName', @(n) true);
parse(p, varargin{:});

%% Prepare input data:
track_files = CollectTrackFiles(data_folder);
track_files = track_files(cellfun(p.Results.FilterByName, track_files));
[~, ~, frames_raw] = xlsread(fullfile(data_folder, 'frames.xlsx'), 'Comments');

track_stats_header = {'Filename', 'Strain', 'Odorant concentration', 'Date',...
    'Total tracks', 'Filtered tracks', 'Mean length', 'Median length',...
    'Total steps', 'Long runs', 'Pirrouetes', 'Max. worms', '% in long runs', ...
    '% in pirouttes'};
track_stats_output = cell(length(track_files) + 1, length(track_stats_header));
track_stats_output(1,:) = track_stats_header;

chemotaxis_stats_header = {'Filename', 'Strain', 'Odorant concentration', 'Date', 'Max. worms', ...
    'First arrival', 'First 10%', 'First 20%', 'First 40%', 'First 80%', ...
    'Arrival at 5', 'Arrival at 10', 'Arrival at 15', 'Arrival at 20'...
    };
chemotaxis_stats_output = cell(length(track_files) + 1, length(chemotaxis_stats_header));
chemotaxis_stats_output(1,:) = chemotaxis_stats_header;

combined_data = containers.Map;

%% Setup output flags and other common constants
% Often we want to regenerate only specific figures. As the code will not
% generate figures if they already exist, these flags allow to force figure
% regeneration.
force_bearings_fig = false;
force_polygon_dynamics_fig = false;
force_att_vec_field = false;
force_distances_fig = false;
force_turn_improvement = false;
force_top_10 = false;
force_top_10_segmented = false;
force_weathervane_fig = false;
    
% In our analysis, we want to separate the circle into 20, and we want to 
% be constitent about it:
bin_edges = -180:20:180;
bin_count = length(bin_edges) - 1;

%% Go over the experiment files:
for track_file_ix = 1:length(track_files)
    load(track_files{track_file_ix})
    
    % Initialize the counters for the track stats:
    output_steps = 0;
    output_pirouettes = 0;
    output_pirouette_steps = 0;
    output_kinks = 0;
    output_long_runs = 0;
    output_long_run_steps = 0;

    % Look in the frames table to find the start and end frames which will
    % be used in the tracks to analyze:
    frames_row = frames_raw(...
        cellfun(@(n) contains(tracker.name, num2str(n)), frames_raw(:,1)) & ...
        (cellfun(@ischar, frames_raw(:,2)) | ...
         ~cellfun(@(n) (~ischar(n) && isnan(n)), frames_raw(:,2))), ...
        2:3);
    start_frame = frames_row{1};
    stop_frame = frames_row{2};
    % In the file, I sometimes specified ranges for the final frame. We
    % want to take the first frame in that range:
    if ischar(stop_frame)
        frame_split = strsplit(stop_frame, '-');
        stop_frame = str2double(frame_split{1});
    end
    
    % Set up shorthands, thresholds and metadata for later use:
    fps = movie_features.fps;
    px_per_mm = movie_features.pixels_per_mm;
    min_track_duration = ceil(30 * fps);
    plate_radius = movie_features.plate{2};
    drop_center = movie_features.drop{1};
    drop_radius = movie_features.drop{2};
    strain = comments.Strain;
    odorant_concentration = comments.OdorantConcentration(end-1:end);
    movie_date = comments.Date;
    label = sprintf('%s - %s - %s', strain, odorant_concentration, movie_date); % To be used in headers
    filename_label = sprintf('%s - %s - %s', strain, odorant_concentration, name); % To be used as output filename base
    
    % How long do we want the pre-pirouette window to be?
    prepirouette_length = round(fps * 60);
    
    % Filter the tracks to analyze:
    filtered_tracks = tracks(...
        ... % Filter by minimum duration, to decrease noise:
        cellfun(@(p) p(end,1) - p(1,1), {tracks.filteredPath}) >= min_track_duration & ...
        ... % Use tracks that start beyond a certain point on the plate, to decrease noise due to circling around the center:
        cellfun(@(p) Dist(drop_center, p(1,2:3)), {tracks.filteredPath}) > 1/4*plate_radius & ...
        ... % Filter by start and end frames:
        cellfun(@(p) p(end,1) > start_frame, {tracks.filteredPath}) & ...
        cellfun(@(p) p(1,1) < stop_frame, {tracks.filteredPath})...
        );
    
    % Sort by track length, mostly for convenience later on:
    track_lengths = cellfun(@sum, {filtered_tracks.filteredStepSizes});
    [sorted_track_lengths, sorted_length_ixes] = sort(track_lengths);
    sorted_tracks = filtered_tracks(sorted_length_ixes);

    % Set up structures for collecting data about tracks:
    pirouette_bearings = [];
    long_run_bearings = [];
    long_run_speeds = [];
    
    % The prepirouette structure collects information about the state of
    % the animal during N seconds before a pirouette - its bearing and
    % speed. 
    % It has a 3D structure with the following indexes:
    % x - prepirouette #
    % y - step ix
    % z - bearing (1) and speed (2)
    prepirouettes = zeros(0, prepirouette_length, 2);
    
    for track_ix = 1:length(sorted_tracks)
        track = sorted_tracks(track_ix);
        fp = track.filteredPath(:,2:3);
        segmented_path = SegmentTrackPath(track, fps);
        
        long_run_ixes = bitand(segmented_path, SegmentMasks.LongRun) > 0;
        long_run_bearings = [long_run_bearings track.bearings(long_run_ixes)'];
        long_run_speeds = [long_run_speeds track.filteredStepSizes(long_run_ixes)' / px_per_mm * fps];
        
        [pirouettes, num_of_pirouettes] = SegmentVector(bitand(segmented_path, SegmentMasks.Pirouette) > 0);
        
        output_steps = output_steps + length(track.filteredStepSizes);
        output_pirouettes = output_pirouettes + sum(pirouettes(:,3) > 0);
        output_pirouette_steps = output_pirouette_steps + sum(bitand(segmented_path, SegmentMasks.Pirouette) > 0);
        output_kinks = output_kinks + sum(bitand(segmented_path, SegmentMasks.Kink) > 0);
        output_long_runs = output_long_runs + sum(pirouettes(:,3) == 0);
        output_long_run_steps = output_long_run_steps + sum(cellfun(@(p) (p(3) == 0) * (p(2) - p(1) + 1), mat2cell(pirouettes, ones(num_of_pirouettes, 1))));
        
        % Analyze the pirouettes:        
        % Start from ix=2 because we must know the pre-pirouette bearing
        % (and if a pirouette is at ix=1 by definition we don't know it).
        for segment_ix = 2:num_of_pirouettes
            if pirouettes(segment_ix,3) == 0
                % Ignore long runs
                continue
            end
            
            pirouette_start = pirouettes(segment_ix,1);
            pirouette_end = pirouettes(segment_ix,2);
            
            long_run_start = pirouettes(segment_ix-1,1);
            long_run_end = pirouettes(segment_ix-1,2);
            if long_run_end - long_run_start >= prepirouette_length
                prepirouette(:,1) = track.bearings(long_run_end-prepirouette_length+1:long_run_end);
                prepirouette(:,2) = track.filteredStepSizes(long_run_end-prepirouette_length+1:long_run_end) / px_per_mm * fps;
                prepirouettes(end+1,:,:) = prepirouette;
            end
            
            % We take the bearing before and after the pirouette at some
            % distance from the pirouette start and end:
            points_count = 6;
            before_origin = fp(pirouette_start - points_count,:);
            bearing_before_pir = GetAngle(drop_center - before_origin, fp(pirouette_start,:) - before_origin);
            
            if segment_ix < num_of_pirouettes
                after_origin = fp(pirouette_end+1,:);
                bearing_after_pir = GetAngle(drop_center - after_origin, fp(pirouette_end+1+points_count,:) - after_origin);
            else
                bearing_after_pir = nan;
            end
            
            pirouette_bearings = [pirouette_bearings; [bearing_before_pir bearing_after_pir]];
        end
    end
    
    % Not all pirouette bearings have an exit bearing, but sometimes we do
    % want to compare the after vs. before bearings:
    matched_bearings = pirouette_bearings(~isnan(pirouette_bearings(:,2)),:);
    
    % For every track, compute the the distances to the drop center
    % per each frame. 
    [distances_per_frame, first_frame_with_track, last_frame_with_track] = GetDistancesPerFrame(tracks, drop_center);
    tracked_worms_per_frame = sum(~isnan(distances_per_frame), 1);
    
    % Plot the absolute speeds against the bearing.
    abs_speed_fig = [fullfile(out_folder, filename_label) '.speeds_vs_bearing.fig'];
    if ~exist(abs_speed_fig, 'file')
        binned_bearings = discretize(long_run_bearings, bin_edges);
        binned_speeds = zeros(bin_count, 1);
        for bin_ix = 1:length(bin_edges) - 1
            binned_speeds(bin_ix) = mean(long_run_speeds(binned_bearings == bin_ix));
        end
        
        figure;
        bar(bin_edges(1:end-1), binned_speeds);
        title(label);
        savefig(abs_speed_fig);
        close;
    end
    
    % Plot the bearing and pirouette info
    fig_filename = [fullfile(out_folder, filename_label) '.bearings_and_pirouettes.fig'];
    if force_bearings_fig || ~exist(fig_filename, 'file')
        normed_bearings = ComputeConditionalTurnProbabilities(long_run_bearings, pirouette_bearings(:,1), bin_edges);
        
        fig_rows = 4;
        fig_cols = 2;
        
        hists = figure;
        
        subplot(fig_rows, fig_cols, 1);
        histogram(long_run_bearings, bin_edges);
        title('# steps vs. B');

        subplot(fig_rows, fig_cols, 2);
        histogram(pirouette_bearings(:,1), bin_edges);
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
        abs_delta_bearings = zeros(bin_count, 1);
        for bin_ix = 1:length(bin_edges) - 1
            abs_delta_bearings(bin_ix) = mean(abs(deltas(binned_bearings == bin_ix)));
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
        xlim([-prepirouette_length/fps 0]);
        title('Prepirouette B vs. t (sec)');
        
        subplot(fig_rows, fig_cols, 8);
        plot(ts, mean(prepirouettes(:,:,2)));
        xlim([-prepirouette_length/fps 0]);
        title('Prepirouette speed vs. t (sec)');

        mtit(hists, label, 'xoff', 0, 'yoff', 0.03, 'zoff', 0);
        savefig(hists, fig_filename);
        close;
    end

    % Compute various statistics about arrival and accumulation near the
    % odorant. We don't plot the polygon dynamics original plot, but draw
    % our own plot, normalized to the maximum number of worms observed.
    worms_in_region = polygonDynamics(tracker, false,...
        'Polygon', circleToPolygon([drop_center plate_radius / 5], 8),...
        'Tracks', tracks);

    normed_worms_in_region = worms_in_region(start_frame:end) / max(tracked_worms_per_frame);
    ts = linspace(0, length(normed_worms_in_region)/fps, length(normed_worms_in_region)) / 60;
    stop_time = (stop_frame-start_frame) / fps / 60;

    first_worm_time = find(normed_worms_in_region > 0) / fps / 60;
    first_worm_time = first_worm_time(1);

    pol_dyn_fig_name = fullfile(out_folder, [filename_label '.norm_polygon_dynamics.fig']);
    if force_polygon_dynamics_fig || ~exist(pol_dyn_fig_name, 'file')
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
    end
    
    % Plot the attraction vector field - where are the worms going to?
    att_vec_field_fig_name = fullfile(out_folder, [filename_label '.attraction_vector_field.fig']);
    if force_att_vec_field || ~exist(att_vec_field_fig_name, 'file') 
        attractionVectorField(tracker, 0, 'Tracks', sorted_tracks');
        title([label ' - ' name]);
        savefig(att_vec_field_fig_name);
        close;
    end
    
    % Draw a plot of average distances of all worms vs. time.
    med_dist_plot_name = fullfile(out_folder, [filename_label '.all_distances.fig']);
    if force_distances_fig || ~exist(med_dist_plot_name, 'file')
        tracked_worms = sum(~isnan(distances_per_frame), 1);
        
        figure;
        hold on;
        title(label);
        mean_dists = mean(distances_per_frame, 1, 'omitnan') / plate_radius;
        median_dists = median(distances_per_frame, 1, 'omitnan') / plate_radius;
        std_dev = std(distances_per_frame, 0, 1, 'omitnan') / plate_radius;
        tracked_worms = tracked_worms / max(tracked_worms); % Normalize
        xs = first_frame_with_track:last_frame_with_track;
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
    
    all_weathervanes = [];
    for track_ix = 1:length(sorted_tracks)
        weathervane = Weathervane(sorted_tracks(track_ix), movie_features);
        if isempty(weathervane)
            continue
        end
        
        weathervane = weathervane( ...
            weathervane(:,3) - start_frame < fps * 60 * 15 & ...
            weathervane(:,4) <= 2 * plate_radius / 3 & ...
            weathervane(:,4) >= 1 * plate_radius / 3, ...
            : ...
        );
        all_weathervanes = [all_weathervanes; weathervane];
    end
    
    turn_improvement_filename = fullfile(out_folder, [filename_label '.turn_improvement.fig']);
    if force_turn_improvement || ~exist(turn_improvement_filename, 'file')
        filtered_bearings = matched_bearings(abs(matched_bearings(:,1)) >= 90,:);
        theta_before = deg2rad(filtered_bearings(:,1));
        theta_after = deg2rad(filtered_bearings(:,2));
        
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
    
    % Plot the absolute speeds of the worms
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
        title({label; 'Speed in all steps'});

        subplot(2, 1, 2);
        histogram(long_run_speeds, edges);
        title('Only long runs');

        savefig(abs_speeds_filename);
        close;
    end
    
    % Store stats:
    accumulation_stats = {};
    for accumulation_threshold = [0.1 0.2 0.4 0.8]
        arrival_time = nan;
        arrival_ix = find(normed_worms_in_region >= accumulation_threshold, 1);
        if ~isempty(arrival_ix)
            arrival_time = datestr(duration(0, ts(arrival_ix), 0), 'MM:SS');
        end
        accumulation_stats{end+1} = arrival_time;
    end

    arrived_percent = [];
    for arrival_time = [5 10 15 20]
        arrived_percent(end+1) = normed_worms_in_region(round(arrival_time * 60 * fps));
    end
    
    chemotaxis_stats_output(track_file_ix+1,:) = {tracker.name, comments.Strain, comments.OdorantConcentration, ...
        movie_date, max(tracked_worms_per_frame), datestr(duration(0, first_worm_time, 0), 'MM:SS'), ...
        accumulation_stats{1}, accumulation_stats{2}, accumulation_stats{3}, accumulation_stats{4}, ...
        arrived_percent(1), arrived_percent(2), arrived_percent(3), arrived_percent(4) ...
    };

    track_stats_output(track_file_ix+1,:) = {tracker.name, comments.Strain, ...
        comments.OdorantConcentration, movie_date,...
        length(tracks), length(sorted_tracks), mean(sorted_track_lengths), median(sorted_track_lengths),...
        output_steps, output_long_runs, output_pirouettes,...
        max(tracked_worms_per_frame), output_long_run_steps / output_steps,...
        output_pirouette_steps / output_steps ...
    };
    
    % Collect data per condition - strain and concentration.
    condition_key = sprintf('%s - %s', comments.Strain, comments.OdorantConcentration);
    if ~isKey(combined_data, condition_key)
        combined_data(condition_key) = struct('all_bearings', [], 'turn_bearings', [], ...
            'weathervane', [], 'prepirouettes', []);
    end
    
    prev_data = combined_data(condition_key);
    combined_data(condition_key) = struct(...
        'all_bearings', [prev_data.all_bearings long_run_bearings], ...
        'turn_bearings', [prev_data.turn_bearings; pirouette_bearings], ...
        'weathervane', [prev_data.weathervane; all_weathervanes], ...
        'prepirouettes', cat(1, prev_data.prepirouettes, prepirouettes) ...
    );
    
    % Print an indicator of progress:
    fprintf('Finished %02d/%02d - %s\n', track_file_ix, length(track_files), label);
end

%% Write out the track and chemotaxis stats:
xlswrite(fullfile(out_folder, 'movie_statistics.xlsx'), track_stats_output);
xlswrite(fullfile(out_folder, 'chemotaxis_statistics.xlsx'), chemotaxis_stats_output);

%% Plot the pooled data per condtiion:
for condition_cell = keys(combined_data)
    condition = condition_cell{1};
    data = combined_data(condition);
    
    % Plot weathervane data:
    all_weathervanes = data.weathervane;
    w_bearing_bins = discretize(all_weathervanes(:,1), bin_edges);
    w_avg_per_bin = zeros(length(bin_edges)-1, 1);
    for bin_ix = 1:length(w_avg_per_bin)
        w_avg_per_bin(bin_ix) = mean(all_weathervanes(w_bearing_bins == bin_ix,2));
    end

    figure;
    bar(bin_edges(1:end-1), w_avg_per_bin)
    title(condition);

    savefig(fullfile(out_folder, [condition '.combined.weathervane.fig']));
    close;

    % Plot turn improvement:
    matched_bearings = data.turn_bearings(~isnan(data.turn_bearings(:,2)),:);
    filtered_bearings = matched_bearings(abs(matched_bearings(:,1)) >= 90 ,:);
    theta_before = deg2rad(filtered_bearings(:,1));
    theta_after = deg2rad(filtered_bearings(:,2));

    polarhistogram(theta_before, 18);
    hold on;
    polarhistogram(theta_after, 18);
    set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');
    title(condition);
    legend('Before', 'After');
    savefig(fullfile(out_folder, [condition '.combined.turn_improvement.fig']));
    close;
    
    %% Plot bearing and pirouette data:
    conditional_probs = ComputeConditionalTurnProbabilities(data.all_bearings, data.turn_bearings(:,1), bin_edges);
    
    fig_rows = 4;
    fig_cols = 2;

    hists = figure;

    subplot(fig_rows, fig_cols, 1);
    histogram(data.all_bearings, bin_edges);
    title('# steps vs. B');

    subplot(fig_rows, fig_cols, 2);
    histogram(data.turn_bearings(:,1), bin_edges);
    title('# pir vs. B_{before}');

    subplot(fig_rows, fig_cols, 3);
    bar(bin_edges(1:end-1), conditional_probs);
    title('P(pir|B_{before})');

    subplot(fig_rows, fig_cols, 4);
    histogram(matched_bearings(:,2), bin_edges);
    title('# pir vs. B_{after}');

    binned_bearings = discretize(matched_bearings(:,1), bin_edges);
    deltas = arrayfun(@(b, a) GetAngle([cosd(b) sind(b)], [cosd(a) sind(a)]), ...
        matched_bearings(:,1), matched_bearings(:,2));
    abs_delta_bearings = zeros(bin_count, 1);
    for bin_ix = 1:bin_count
        abs_delta_bearings(bin_ix) = mean(abs(deltas(binned_bearings == bin_ix)));
    end

    subplot(fig_rows, fig_cols, 5);
    histogram(deltas, bin_edges);
    title('# pir vs. dB');

    subplot(fig_rows, fig_cols, 6);
    bar(bin_edges(1:end-1), abs_delta_bearings);
    title('Abs mean dB vs. B_{before}');

    [~, prepirouette_length, ~] = size(data.prepirouettes);
    ts = linspace(-prepirouette_length/fps, 0, prepirouette_length);

    subplot(fig_rows, fig_cols, 7);
    plot(ts, mean(abs(data.prepirouettes(:,:,1))));
    hold on;
    plot(ts, std(abs(data.prepirouettes(:,:,1))));
    xlim([-prepirouette_length/fps 0]);
    title('Prepirouette B vs. t (sec)');

    subplot(fig_rows, fig_cols, 8);
    plot(ts, mean(data.prepirouettes(:,:,2)));
    xlim([-prepirouette_length/fps 0]);
    title('Prepirouette speed vs. t (sec)');

    mtit(hists, condition, 'xoff', 0, 'yoff', 0.03, 'zoff', 0);
    savefig(hists, fullfile(out_folder, [condition '.combined.bearings_and_pirouettes.fig']));
    close;
end

end

function [mat_files] = CollectTrackFiles(folder)

mat_files = {};

listing = dir(folder);

for ix = 1:numel(listing)
    entry = listing(ix);
    if entry.isdir && ~strcmp(entry.name, '.') && ~strcmp(entry.name, '..')
        mat_files = horzcat(mat_files, CollectTrackFiles(fullfile(entry.folder, entry.name)));
    elseif contains(entry.name, '.preprocessed_tracks.mat')
        mat_files{end+1} = fullfile(entry.folder, entry.name);
    end
end

end

function conditional_bins = ComputeConditionalTurnProbabilities(all_bearings, bearings_before_turn, bin_edges)

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
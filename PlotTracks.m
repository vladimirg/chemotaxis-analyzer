function out_fig = PlotTracks(tracks, movie_features, varargin)

p = inputParser;
addParameter(p, 'IndexLabels', -1);
addParameter(p, 'MarkIndexes', {});
addParameter(p, 'MarkSymbol', 'x');
addParameter(p, 'SegmentWithFps', -1);
addParameter(p, 'OmitIndexing', true);
addParameter(p, 'OmitKinks', true);
addParameter(p, 'OmitSharpTurns', true);
addParameter(p, 'Smoothing', 'None'); % Can also take 'Add' and 'Only'
%addParameter(p, 'TimeSeries', 3); % In FPS
parse(p, varargin{:});

out_fig = figure;
hold on;
viscircles([movie_features.plate{1}; movie_features.drop{1}], [movie_features.plate{2}; movie_features.drop{2}]);
for ix = 1:length(tracks)
    track = tracks(ix);
    if isfield(track, 'filteredPath')
        fp = track.filteredPath(:,2:3);
    else
        fp = track.path(:,2:3);
    end
    
    if strcmp(p.Results.Smoothing, 'Only')
        plot(smooth(fp(:,1)), smooth(fp(:,2)));
        continue
    end
    
    plot(fp(:,1), fp(:,2));

    if p.Results.IndexLabels ~= -1
        arrayfun(@(i) text(fp(i, 1), fp(i, 2), num2str(i)), 1:p.Results.IndexLabels:length(fp));
    end
    
    if ~isempty(p.Results.MarkIndexes)
        if ~iscell(p.Results.MarkIndexes)
            mark_indexes = p.Results.MarkIndexes;
        else
            mark_indexes = p.Results.MarkIndexes{ix};
        end
        
        plot(fp(mark_indexes,1), fp(mark_indexes,2), p.Results.MarkSymbol);
    end
    
    if p.Results.SegmentWithFps > 0
        segments_vector = SegmentTrackPath(track, p.Results.SegmentWithFps);
        
        % Paint long runs green
        [long_runs, long_run_count] = SegmentVector(bitand(segments_vector, SegmentMasks.LongRun) > 0);
        for segment_ix = 1:long_run_count
            if ~long_runs(segment_ix,3)
                continue
            end
            
            start_point_ix = long_runs(segment_ix,1);
            end_point_ix = long_runs(segment_ix,2) + 1;
            
            plot(fp(start_point_ix:end_point_ix,1), ...
                fp(start_point_ix:end_point_ix,2), ...
                'g');
            
            start_point = fp(start_point_ix,:);
            end_point = fp(end_point_ix,:);
            
            if ~p.Results.OmitIndexing
                text(start_point(1), start_point(2), num2str(start_point_ix));
                text(end_point(1), end_point(2), num2str(end_point_ix));
            end
        end
        
        % Paint pirouttes red
        % Or, anything that's between a long run!
%         [pirouttes, pirouette_count] = SegmentVector(bitand(segments_vector, double(bitcmp(uint8(SegmentMasks.LongRun)))) > 0);
        for segment_ix = 1:long_run_count
            if long_runs(segment_ix,3)
                continue
            end
            
            start_point_ix = long_runs(segment_ix,1);
            end_point_ix = long_runs(segment_ix,2) + 1;
            
            plot(fp(start_point_ix:end_point_ix,1), ...
                fp(start_point_ix:end_point_ix,2), ...
                'r');
            
            pirouette_start = start_point_ix;
            pirouette_end = end_point_ix - 1;
            
            points_count = 6;
            
            if segment_ix > 1
                before_origin = fp(pirouette_start - points_count,:);
                plot([before_origin(1) fp(pirouette_start,1)],...
                    [before_origin(2) fp(pirouette_start,2)],...
                    'b');
            end
            
            if segment_ix < long_run_count
                after_origin = fp(pirouette_end+1,:);
                plot([after_origin(1) fp(pirouette_end+1+points_count,1)],...
                    [after_origin(2) fp(pirouette_end+1+points_count,2)],...
                    'b');
            end
        end
        
        if ~p.Results.OmitSharpTurns
            % Mark sharp turns as Xs
            sharp_turns = bitand(segments_vector, SegmentMasks.SharpTurn) > 0;
            plot(fp(sharp_turns,1), fp(sharp_turns,2), 'x');
        end
        
        if ~p.Results.OmitKinks
            % Mark kinks as Os
            kinks = bitand(segments_vector, SegmentMasks.Kink) > 0;
            plot(fp(kinks,1), fp(kinks,2), 'o');
        end
    end

    if p.Results.SegmentWithFps < 0 || p.Results.OmitIndexing
        % Mark start and end of tracks:
        text(fp(1,1), fp(1,2), '1');
        [steps, ~] = size(fp);
        text(fp(end,1), fp(end,2), num2str(steps));
    end
    
    if strcmp(p.Results.Smoothing, 'Add')
        plot(smooth(fp(:,1)), smooth(fp(:,2)));
    end
    
    pbaspect([1 1 1]);
end

end
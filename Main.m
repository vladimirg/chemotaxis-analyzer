% FixVideoPath('F:\VladMovies\Tracked movies\Tracks', 'F:\VladMovies\Tracked movies');
% FeatureExtractor('F:\VladMovies\Tracked movies\Tracks', 'F:\tmp\test tracks');
% TrackAnalyzer('F:\VladMovies\Tracked movies\Tracks', 'F:\tmp\test tracks\MovieFeatures.mat', 'F:\tmp\test tracks\');
% AnalyzePreprocessedTracks('F:\tmp\test tracks', 'F:\tmp\test tracks\figs-feb-tsevi', ...
%     'FilterByName', @(n) contains(n, 'Feb')...
%     );

AnalyzePreprocessedTracks('F:\tmp\test tracks\tsevi test', 'F:\tmp\test tracks\tsevi test - out');


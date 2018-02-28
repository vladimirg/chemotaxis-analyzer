function [mat_files] = GetTrackingInfoFiles(folder)

mat_files = {};

listing = dir(folder);

for ix = 1:numel(listing)
    entry = listing(ix);
    if entry.isdir && ~strcmp(entry.name, '.') && ~strcmp(entry.name, '..')
        mat_files = horzcat(mat_files, GetTrackingInfoFiles(fullfile(entry.folder, entry.name)));
    elseif strcmp(entry.name, 'TrackingInfo.mat')
        mat_files{end+1} = fullfile(entry.folder, entry.name);
    end
end

end
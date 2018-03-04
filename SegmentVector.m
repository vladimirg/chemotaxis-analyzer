function [segments, varargout] = SegmentVector(vector)
% Returns the segments, and optionally the number of segments (in the
% second parameter). Segments are defined as consecutive runs of the same
% value. The output is Nx3, where N is the number of segments, the first
% column is the first index of the segment, the second column is the
% last index of the segment, and the third column is the value of the
% index. For example, [0 1 1 0] will be segmented into:
% [1 1 0]
% [2 3 1]
% [4 4 0]

if isempty(vector)
    segments = [];
    return
end

segments = zeros(sum(vector(2:end) ~= vector(1:end-1)) + 1, 3);
prev_state = vector(1);
prev_ix = 1;
segment_ix = 1;
for ix = 1:length(vector)
    state = vector(ix);
    if state ~= prev_state
        segments(segment_ix,:) = [prev_ix, ix-1, prev_state];
        prev_state = state;
        segment_ix = segment_ix + 1;
        prev_ix = ix;
    end
end

segments(end,:) = [prev_ix, length(vector), state];

if nargout == 2
    [num_of_segments, ~] = size(segments);
    varargout{1} = num_of_segments;
end

end
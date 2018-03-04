function turning_rate = GetTurningRate(track, fps)
% Get the turning rate (first derivate) of the track based on the absolute
% bearings. Absolute bearings are relative to a fixed vector (say, [1 0]),
% and not the odorant center (which is a fixed point).

% TODO: do we really need to smooth the bearings now that we smooth the
% path?

if isfield(track, 'filteredPath')
    path = track.filteredPath;
else
    path = track.path;
end

ts = (path(2:end,1) - path(2,1)) * (1/fps);
fixed_bearings = CircularizeBearings(track.absoluteBearings);
smoothed_bearings = smooth(fixed_bearings);
turning_rate = diffxy(ts, smoothed_bearings);

% This code allows visualization of the relevant data, and is useful for
% debugging.
% xs = 1:length(turning_rate);
% figure;
% hold on;
% plot(xs, ones(length(smoothed_bearings),1) * 50, '--k')
% plot(xs, ones(length(smoothed_bearings),1) * -50, '--k')
% 
% % plot(ts, track.absoluteBearings, 'g');
% % plot(ts, turning_rate, 'r');
% plot(xs, smoothed_bearings);
% plot(xs, turning_rate);


end
function [] = viewTrack( track )
    track.vt.viewTracks(track, [track.path(1,1) track.path(end,1)]);
end


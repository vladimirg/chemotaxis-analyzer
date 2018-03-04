classdef SegmentMasks
    properties (Constant)
        ShortRun = 1;
        LongRun = 2;
        SharpTurn = 4;
        Kink = 8; % TODO: not used anymore; consider removing?
        Pirouette = 16;
        Reversal = 32; % TODO: not used right now.
        OmegaTurn = 64; % TODO: not used right now.
    end
end
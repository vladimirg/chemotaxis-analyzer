classdef SegmentMasks
    properties (Constant)
        ShortRun = 1;
        LongRun = 2;
        SharpTurn = 4;
        Kink = 8;
        Pirouette = 16; % Two sharp turns separated by at least 2 seconds
        Reversal = 32; % Two sharp turns under 2 seconds in total, without change of bearing
        OmegaTurn = 64; % A single sharp turn, but not a reversal - <= 160 degrees.
    end
end
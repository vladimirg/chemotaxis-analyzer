function angle = GetAngle(v1, v2)
% Let v1 = [x1,y1] and v2 = [x2,y2], this function returns the angle in
% degrees between the vectors as measured in a counterclockwise direction
% from v1 to v2. If that angle would exceed 180 degrees, then the angle is
% measured in the clockwise direction but given a negative value.
%
% For computing absolute angles in a unit circle where r=[1 0] is the
% reference vecotr, set v1=r. So, GetAngle(r, [0.5 0.5]) = 45, and
% GetAngle(r, [0.5 -0.5]) = -45.
%
% When dealing with relative angles, remember that the clockwise direction
% from v1 to v2 is negative, while the counterclockwise is positive. As a
% corollary, a positive angle means that v2 is the right of v1, while a
% negative angle means that v2 is to the left of v1 (always relative to v1,
% of course).

% Code borrowed from:
% https://www.mathworks.com/matlabcentral/answers/180131-how-can-i-find-the-angle-between-two-vectors-including-directional-information

x1 = v1(1); y1 = v1(2);
x2 = v2(1); y2 = v2(2);

angle = atan2d(x1*y2-y1*x2,x1*x2+y1*y2);

end
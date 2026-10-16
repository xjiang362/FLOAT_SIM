function zCorners = cornerHeightSensor(x, params)
%CORNERHEIGHTSENSOR Simulated 4 corner height (world z) measurements --
%what a real gap/height sensor at each corner would read. Uses the full
%rotation (not a small-angle approximation), since this represents the
%physical measurement, not the controller's estimate of it.

r = x(1:3);
q = x(7:10);
R = quat2rotmCustom(q);

cornersBody = [params.mixer.rx'; params.mixer.ry'; zeros(1, 4)];   % 3x4
cornersWorld = R * cornersBody + r;                                % 3x4

zCorners = cornersWorld(3, :)';   % 4x1

end

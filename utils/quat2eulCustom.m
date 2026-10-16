function eul = quat2eulCustom(q)
%QUAT2EULCUSTOM ZYX (yaw-pitch-roll) Euler angles [rad] from a quaternion
%[w;x;y;z]. Used only for readable plots, not for propagation.

q = quatNormalize(q);
w = q(1); x = q(2); y = q(3); z = q(4);

roll  = atan2(2*(w*x + y*z), 1 - 2*(x^2 + y^2));
pitch = asin(max(-1, min(1, 2*(w*y - z*x))));
yaw   = atan2(2*(w*z + x*y), 1 - 2*(y^2 + z^2));

eul = [yaw; pitch; roll];

end

function q = eul2quatCustom(eul)
%EUL2QUATCUSTOM Quaternion [w;x;y;z] from ZYX (yaw-pitch-roll) Euler
%angles [rad], eul = [yaw; pitch; roll].

yaw = eul(1); pitch = eul(2); roll = eul(3);

cy = cos(yaw/2);   sy = sin(yaw/2);
cp = cos(pitch/2); sp = sin(pitch/2);
cr = cos(roll/2);  sr = sin(roll/2);

q = [cr*cp*cy + sr*sp*sy;
     sr*cp*cy - cr*sp*sy;
     cr*sp*cy + sr*cp*sy;
     cr*cp*sy - sr*sp*cy];

q = quatNormalize(q);

end

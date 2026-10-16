function plotTrapResults(log, setpoint)
%PLOTTRAPRESULTS Diagnostic plots for a trapDeviceEOM simulation log:
%position (with trap center at the origin) / altitude tracking, attitude,
%thruster commands, and trap stiffness commands.

t = log.t;
r = log.x(1:3, :);
q = log.x(7:10, :);
w = log.x(11:13, :);

N = numel(t);
eul = zeros(3, N);
for k = 1:N
    eul(:, k) = quat2eulCustom(q(:, k));
end
eulDeg = rad2deg(eul);

figure('Name', 'Trap Device - Position', 'Color', 'w');
subplot(2,1,1);
plot(t, r(1,:), t, r(2,:), 'LineWidth', 1.4); hold on;
yline(0, 'k:');
grid on; ylabel('Position [m]'); legend('x','y','trap center','Location','best');
title('x, y (passively trapped about the origin)');

subplot(2,1,2);
plot(t, r(3,:), 'LineWidth', 1.4); hold on;
yline(setpoint.z, 'r--');
grid on; xlabel('Time [s]'); ylabel('Altitude z [m]');
legend('z','z_{desired}','Location','best');
title('Altitude (actively controlled)');

figure('Name', 'Trap Device - Attitude', 'Color', 'w');
subplot(2,1,1);
plot(t, eulDeg, 'LineWidth', 1.4);
grid on; ylabel('Angle [deg]');
legend('yaw (passive)','pitch (active)','roll (active)','Location','best');
title('Attitude (ZYX Euler)');

subplot(2,1,2);
plot(t, rad2deg(w), 'LineWidth', 1.4);
grid on; xlabel('Time [s]'); ylabel('Body rate [deg/s]');
legend('wx','wy','wz','Location','best');
title('Angular Velocity');

figure('Name', 'Trap Device - Actuators', 'Color', 'w');
subplot(2,1,1);
plot(t, log.T, 'LineWidth', 1.2);
grid on; ylabel('Thrust [N]');
legend('T1','T2','T3','T4','Location','best');
title('Corner Thruster Commands');

subplot(2,1,2);
plot(t, log.u, 'LineWidth', 1.2);
grid on; xlabel('Time [s]'); ylabel('Trap gain u [-]');
legend('u_x','u_y','u_{yaw}','Location','best');
title('Trap Stiffness Commands (bang-bang schedule)');

end

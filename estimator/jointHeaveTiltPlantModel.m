function [A, B, C, stateNames] = jointHeaveTiltPlantModel(params)
%JOINTHEAVETILTPLANTMODEL Low-fidelity LINEAR joint model of the heave/
%roll/pitch subsystem (x, y, yaw are not part of this model -- they stay
%with the separate, nonlinear trap mechanism).
%
%   State x = [z; vz; roll; wx; pitch; wy]        (6x1)
%   Input u = [Fz_pert; taux; tauy]                (3x1)
%             Fz_pert = Fz - m*g, the PERTURBATION away from hover thrust
%             (matches heaveController.m's convention: the regulator
%             commands this perturbation, gravity feedforward is added
%             back separately) -- NOT the raw commanded Fz.
%   Output y = 4 corner height readings           (4x1)
%
%   Continuous-time: xdot = A*x + B*u,  y = C*x
%
%   A, B are BLOCK-DIAGONAL -- heave, roll, pitch are three independent
%   double integrators with NO cross-coupling terms. This is not an
%   approximation of convenience; it is exact to first order around the
%   hover/level equilibrium, because the only physical coupling
%   mechanism (the gyroscopic term w x (Iw) in Euler's equation) is
%   quadratic in the rate perturbations -- every term is a product of
%   two different rate components (e.g. wy*wz) -- and both the value and
%   every first partial derivative of such a product vanish at w=0. A
%   linearization around hover cannot recover coupling that is identically
%   zero at the point being linearized around.
%
%   C is NOT block-diagonal: each of the 4 corner sensors reads a linear
%   combination of z, roll, AND pitch simultaneously,
%       z_i = z + ry_i*roll - rx_i*pitch,
%   through plain rigid-body kinematics (same geometry as
%   trapDeviceParams.m's cornerSensor.M, just interleaved with zero
%   columns for the unobserved rate states). This is the coupling that
%   actually matters once sensor data arrives serially/asynchronously
%   with per-corner delay: the plant never mixes the three channels, but
%   a single corner reading always constrains a combination of all three.

m   = params.m;
Ixx = params.I(1, 1);
Iyy = params.I(2, 2);
rx  = params.mixer.rx;
ry  = params.mixer.ry;

blk = [0, 1; 0, 0];
A = blkdiag(blk, blk, blk);

B = zeros(6, 3);
B(2, 1) = 1 / m;     % Fz_pert -> vzdot
B(4, 2) = 1 / Ixx;   % taux    -> wxdot
B(6, 3) = 1 / Iyy;   % tauy    -> wydot

C = zeros(4, 6);
C(:, 1) = 1;     % every corner sees z directly
C(:, 3) = ry;    % roll contributes ry_i*roll to corner i
C(:, 5) = -rx;   % pitch contributes -rx_i*pitch to corner i

stateNames = {'z', 'vz', 'roll', 'wx', 'pitch', 'wy'};

end

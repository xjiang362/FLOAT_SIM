function log = simulateTrapLevitation(plant, K, x0, tEnd, varargin)
%SIMULATETRAPLEVITATION Closed-loop levitation simulation.
%
%   log = simulateTrapLevitation(plant, K, x0, tEnd)
%   log = simulateTrapLevitation(plant, K, x0, tEnd, 'delay', false)
%
% Integrates the full 12-state nonlinear plant with forward Euler at
% plant.sim.dtEuler, running the discrete controller every plant ISR tick.
% The loop is capacitance -> sensing -> height -> estimator -> regulator
% -> force -> voltage, so every stage built so far is exercised.
%
% INPUTS
%   plant  struct from trapDeviceGeometry()
%   K      controller from trapControlDesign()
%   x0     12-by-1 [pose(6); rate(6)], metres and radians
%   tEnd   [s] simulation length
%
% NAME/VALUE
%   'delay'   logical, include the acquisition delay (default true)
%   'ref'     3-by-1 setpoint [rz; thx; thy] (default zeros)
%   'Vtia'    [V] fixed bias on the TIA pair (default 40)
%   'trueState' logical, bypass the estimator and feed the regulator the
%            exact plant state (default false). Not implementable on
%            hardware -- it is a diagnostic that removes the estimator from
%            the loop so the zero-order hold can be studied on its own.
%
% OUTPUT log
%   .t .x .u .f .zc .xhat .margin .contact .Ceff
%
% The four sensor channels arrive with different ages when 'delay' is on
% (1130, 950, 770, 590 us), so the controller is fed a set of readings
% that refer to four different past poses. With 'delay' off every channel
% is read at the current instant, which is the design assumption.
%
% Integration stops early if the plate contacts a wall; log.contact marks
% the index.

ip = inputParser;
ip.addParameter('delay', true);
ip.addParameter('ref', zeros(3,1));
ip.addParameter('Vtia', 40);
ip.addParameter('Vmax', 80);
ip.addParameter('trueState', false);
ip.parse(varargin{:});
useDelay = ip.Results.delay;
ref = ip.Results.ref(:);
Vtia = ip.Results.Vtia;
Vmax = ip.Results.Vmax;
trueState = ip.Results.trueState;

dt   = plant.sim.dtEuler;
Tisr = K.Tisr;
nSub = round(Tisr/dt);
if abs(nSub*dt - Tisr) > 1e-12
    error('simulateTrapLevitation:step','dtEuler must divide Tisr evenly.');
end
nTick = floor(tEnd/Tisr);

% Seed the acquisition history from the initial pose. Without this the
% sampler cannot serve a delayed sample for the first ceil(maxAge/Tisr)
% ticks, the loop silently runs delay-free over that window, and the delay
% then appears in one step -- a disturbance large enough to destabilise the
% plant on its own. See the note in trapSamplerInit.
s = trapSamplerInit(plant, Tisr, 'seedPose', x0(1:6));

% sensing-point coordinates, for seeding the estimator
sInfoPt = zeros(4,2);
for q = 1:4
    ea = plant.term(plant.idx.sensePair(q,1));
    eb = plant.term(plant.idx.sensePair(q,2));
    sInfoPt(q,:) = [(ea.Xc+eb.Xc)/2, (ea.Yc+eb.Yc)/2];
end

x    = x0(:);
xhat = zeros(6,1);                 % estimator state
xi   = zeros(3,1);                 % integrator state
f    = zeros(3,1);
u    = trapWrenchToVoltage(plant, f, zeros(6,1), 'Vtia', Vtia, 'Vmax', Vmax);

N = nTick*nSub + 1;
log.t = zeros(N,1);  log.x = zeros(N,12); log.margin = zeros(N,1);
log.u = zeros(nTick,6); log.f = zeros(nTick,3);
log.zc = nan(nTick,4); log.xhat = zeros(nTick,6); log.Ceff = nan(nTick,4);
log.tTick = zeros(nTick,1);
log.contact = [];

n = 1; log.t(1)=0; log.x(1,:)=x';
for kTick = 1:nTick
    t = (kTick-1)*Tisr;
    pose = x(1:6);

    % The substep guard below checks the pose BEFORE each Euler step, so a
    % pose that goes invalid on the last substep of a tick is not caught
    % until here. Without this check the sensing call would be the first
    % thing to see it and would throw trapCapacitance:contact instead of
    % reporting contact cleanly.
    [okTick, mTick] = checkPoseValid(plant, pose);
    if ~okTick
        log.contact = n; log.contactTime = t;
        break
    end

    % Everything from here to the end of the tick evaluates the capacitance
    % model, which refuses a non-positive gap. The corner test above is the
    % normal guard, but it is evaluated at the plate corners while
    % trapCapacitance tests each electrode, and trapForceToVoltage evaluates
    % at the clamped ESTIMATE rather than the true pose. A badly tuned
    % regulator can therefore still drive a stage outside the model domain.
    % Leaving the model domain IS contact, so report it as such instead of
    % letting the exception escape the harness.
    try

    % ---------- sensing ----------
    [s, CeDelayed, sInfo] = trapSamplerStep(s, pose, t);
    if useDelay && sInfo.valid
        Ce = CeDelayed;
    else
        Ce = sInfo.CeffNow;                 % no delay, or buffer still filling
    end
    [zc, ~] = trapCapToHeight(plant, Ce);

    ys = zc;                                % 4 heights

    % Seed the estimator from the first reading rather than from zero.
    % Starting at zero makes the first correction a step of the full
    % initial displacement, which the fast estimator turns into a large
    % spurious velocity.
    if kTick == 1
        M0 = [ones(4,1), sInfoPt(:,2), -sInfoPt(:,1)];
        p0 = M0\ys;
        xhat([1 3 5]) = p0;
    end

    % ---------- which state does the controller act on? ----------
    % 'trueState' bypasses the estimator entirely and feeds the regulator
    % the exact plant state. This is not implementable -- rz, thx, thy and
    % all three rates are not measured -- but it separates the two failure
    % mechanisms cleanly. With the true state the pose-estimate lag is zero,
    % so the force-to-voltage inversion cancels the electrostatic stiffness
    % exactly at the instant it is computed. Whatever instability survives
    % is then due only to holding that voltage for one ISR.
    if trueState
        xUse = [x(3); x(9); x(4); x(10); x(5); x(11)];
    else
        xUse = xhat;
    end

    % ---------- regulator with integral action ----------
    yc = K.Cc*xUse;
    f  = -K.Kx*xUse - K.Ki*xi;
    xi = xi + Tisr*(ref - yc);              % integrate the error

    % ---------- force -> voltage, using the LIVE state/estimate ----------
    poseEst = zeros(6,1);
    poseEst(3) = xUse(1); poseEst(4) = xUse(3); poseEst(5) = xUse(5);
    u = trapWrenchToVoltage(plant, f, poseEst, 'Vtia', Vtia, 'Vmax', Vmax);

    log.tTick(kTick)=t; log.u(kTick,:)=u'; log.f(kTick,:)=f';
    log.zc(kTick,:)=zc'; log.xhat(kTick,:)=xUse'; log.Ceff(kTick,:)=Ce';

    % ---------- propagate the plant, holding u ----------
    for j = 1:nSub
        % stop before the pose leaves the model's domain, so the capacitance
        % stage is never asked for a non-positive gap
        [okNow, mNow] = checkPoseValid(plant, x(1:6));
        if ~okNow
            log.contact = n; log.contactTime = t+(j-1)*dt;
            break
        end
        [xd, dg] = trapDynamics(plant, x, u);
        x = x + dt*xd;
        n = n+1;
        log.t(n)=t+j*dt; log.x(n,:)=x'; log.margin(n)=mNow;
    end
    % ---------- estimator, PREDICTOR form ----------
    % L was designed by place(Ad'',Csd'',zE)'', which is the gain for
    %   xhat[k+1] = Ad*xhat[k] + Bd*u[k] + L*(y[k] - Csd*xhat[k])
    % Applying it as a separate correct-then-predict pair double-counts it.
    if ~trueState
        xhat = K.Ad*xhat + K.Bd*f + K.L*(ys - K.Csd*xhat);
    end

    catch ME
        if strcmp(ME.identifier,'trapCapacitance:contact')
            log.contact = n; log.contactTime = t;
            log.contactReason = ME.message;
            break
        else
            rethrow(ME)
        end
    end

    if ~isempty(log.contact), break; end
end

log.t = log.t(1:n); log.x = log.x(1:n,:); log.margin = log.margin(1:n);
log.nTick = kTick;
log.delay = useDelay;
end



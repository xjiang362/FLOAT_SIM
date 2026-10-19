# FLOATS_SIM

Simulation of the FLOATS electrostatic levitation trap: a 5.95 × 5.95 × 0.35 mm
silicon proof mass (28.87 mg) floating between electrode arrays, held by
voltages on 20 electrodes arranged as 10 differential pairs.

The plant is open-loop **unstable**. At fixed voltage the electrostatic
stiffness is negative,

```
dFz/drz = 2W/g0 = +37.76 N/m   ->   p = sqrt(2g/g0) = 1143.7 rad/s = 182 Hz
```

so the proof mass must be actively stabilised in heave and both tilts. The
remaining three degrees of freedom (`rx`, `ry`, `thz`) are passively trapped
and uncontrolled.

Requires **MATLAB** with the **Control System Toolbox** (`c2d`, `ss`, `place`,
`dlqr`, `acker`). Developed on R2025b.

---

## Quick start

```matlab
addpath(genpath(pwd));
verifyTrapPipeline            % end-to-end sanity checks on every stage
plotTrapLevitationDesign      % closed-loop response at the current design point
```

---

## Pipeline structure

Data flows `pose -> capacitance -> sensing -> estimate -> control -> voltage ->
wrench -> pose`. Each stage is one file.

```
config/
  trapDeviceGeometry.m      single source of truth: electrode footprints and
                            centroids, mass properties, gaps, squeeze-film
                            damping, readout timing. Everything else reads
                            geometry from here rather than hard-coding it.

dynamics/
  trapCapacitance.m         C(pose) and dC/dq for all 20 terminals. Pure
                            parallel plate, no fitted constants. Lateral
                            motion enters through overlap truncation at the
                            plate rim, which is what produces lateral force.
  trapWrench.m              energy method at fixed potential,
                            Q_j = 0.5 * sum_i V_i^2 dC_i/dq_j
  trapOpenLoopStiffness.m   fixed-voltage stiffness; with 'afterInversion'
                            returns what survives the voltage inversion
  trapDynamics.m            12-state rigid-body EOM (pose + rates)
  checkPoseValid.m          corner-wise wall clearance

estimator/
  trapSensing.m             Y-delta reduction through the floating proof mass
                            -> four effective capacitances
  trapCapToHeight.m         Ceff -> four corner heights -> (rz, thx, thy)
  trapSamplerInit/Step.m    acquisition model: channels converted one at a
                            time, then shipped over the bus, so each arrives
                            stale by a different amount

control/
  trapControlDesign.m       LQR + integrators on the exact sampled map,
                            observer on the stiffness-retaining map
  trapWrenchToVoltage.m     allocates a commanded wrench [Fz; taux; tauy]
                            across the five drive pairs -> 6 voltages
  trapForceToVoltage.m      older per-pair force inversion, kept for the
                            round-trip verifier

simulateTrapLevitation.m    closed-loop simulation tying all of the above
                            together

presentation/               verification scripts, COMSOL comparison studies
                            and the raw COMSOL exports
```

Current design point (defaults in `trapDeviceGeometry.m` and
`trapControlDesign.m`):

| parameter | value |
|---|---|
| `Tisr` | 500 µs |
| `dtChannel` | 125 µs |
| `tI2C` | 500 µs |
| `dtEuler` | 10 µs |
| `lambdar` | 200 Hz |
| `lambdae` | 500 Hz |

---

## Running the ZOH / feedback-linearisation comparison

This is the central control result in the project, so it is documented in
full.

### What is being compared

The wrench allocation in `trapWrenchToVoltage.m` cancels the electrostatic
stiffness exactly, which makes it tempting to design the regulator on a
double integrator. **That is wrong.** Within one ISR the voltage is frozen,
but the commanded acceleration is not, because

```
u(t) = eps0*A*V[k]^2 / (2m (d0 - rz(t))^2) - g
```

still depends on `rz(t)`, which keeps moving. Expanding about `rz[k]` gives
`u(t) ~ u[k] + p^2 * rz(t)`, so the exact one-tick map is

```
            naive (uncompensated)        compensated (exact sampled map)
    A_d  =  [ 1   T ]                    [ 1          sinh(pT)/p ]
            [ 0   1 ]                    [ 0          cosh(pT)   ]

    B_d  =  [ T^2/2 ]                    [ (cosh(pT)-1)/p^2 ]
            [ T     ]                    [ sinh(pT)/p       ]
```

The inversion cancels the **position** term (the zero in the corner) but the
**velocity** is still amplified by `cosh(pT)` every tick — 1.168 at
`Tisr = 500 µs`, 1.338 at 700 µs. The naive model assumes that factor is 1,
under-damps every correction, and crashes into a wall.

`p = sqrt(2g/g0) = 1143.7 rad/s`, derived in
`dynamics/trapOpenLoopStiffness.m` and confirmed against the 6-DOF model.

### Script 1 — heave only, closest to the derivation

```matlab
addpath(genpath(pwd));
simulateHeaveFeedbackLin
```

**File:** `presentation/simulateHeaveFeedbackLin.m`
**Output:** `presentation/HeaveFeedbackLin_rz.png`

A 1-DOF study at `Tisr = 700 µs`, `lambdar = 100 Hz`, `lambdae = 300 Hz`,
released from `rz = -5 µm`. Other DOF start at zero and stay there by
symmetry. Two subplots: true vs estimated `rz`, and the commanded
acceleration as a ZOH staircase. Expected console output:

```
design model                                 cosh(pT)   peak[um]           fate
naive DT  [1 T; 0 1]                           1.0000    -15.002 contact 8.7 ms
compensated  [1 sinh(pT)/p; 0 cosh(pT)]        1.3379     -5.000        settles
```

To change the operating point, edit near the top of the file:

```matlab
T       = 700e-6;        % ISR period
lambdar = 2*pi*100;      % regulator bandwidth
lambdae = 2*pi*300;      % observer bandwidth
rz0     = -5e-6;         % release height
```

> The observer gain is deliberately placed on the **stiffness-retaining** map
> `[cosh pT, sinh pT/p; p sinh pT, cosh pT]` in both cases, not on each
> regulator's own design model. The inversion cancels the stiffness at the
> *estimated* pose, so the estimation error sees it uncancelled:
> `e[k+1] = (A_s - L*C) e`. Placing `L` on the regulator model makes the
> observer diverge and both designs fail, which would confound the
> comparison.

### Script 2 — full 6-DOF plant

```matlab
addpath(genpath(pwd));
compareDesignModels      % runs both designs, saves data and a 6-panel figure
plotRzZoom               % heave only, first 20 ms, from the saved data
```

**Files:** `presentation/compareDesignModels.m`, `presentation/plotRzZoom.m`
**Outputs:** `CompareDesignModels.png`, `CompareDesignModels_rzZoom.png`,
`CompareDesignModels.mat`

Runs the full nonlinear plant from `(-40, -40, -13) µm` with tilt, at
`Tisr = 500 µs`, `lambdar = 50 Hz`. Expected:

```
  exact sampled map    survived 1000 ms: rz=+6.94e-01 um ...
  double integrator    CONTACT at 7.6 ms
```

`compareDesignModels` saves `CompareDesignModels.mat`, so `plotRzZoom` can be
re-run on the same data without re-simulating. The `.mat` is gitignored as
regenerable.

Run `compareDesignModels` at `lambdar = 200 Hz` and **both** designs fail at
1.1 ms — 200 Hz cannot capture from −13 µm whatever model the gains came
from, so the comparison must be made at a bandwidth where the loop works.
`lamR` is set near the top of the file for this reason.

### Script 3 — adding acquisition delay

```matlab
addpath(genpath(pwd));
compareFeedbackLinRz
```

**File:** `presentation/compareFeedbackLinRz.m`
**Outputs:** `CompareFeedbackLin_rz.png` (3 traces),
`CompareFeedbackLin_rz_noDelay.png` (2 traces)

Same comparison plus a third trace with the acquisition delay enabled:

```
case                                                   peak[um]               fate
no delay, FL NOT compensated (double integrator)         12.859     contact 7.6 ms
no delay, FL compensated (exact sampled map)              7.032    survives 600 ms
DELAY, FL compensated                                     8.811     contact 2.1 ms
```

The two failures have different causes. The uncompensated one is a **model
error**, fixed for free in software. The delayed one is a **latency limit**:
no design model or bandwidth addresses it.

### Switching the design model yourself

`trapControlDesign` takes a `'designModel'` option:

```matlab
p = trapDeviceGeometry();
Kc = trapControlDesign(p, 'designModel', 'exact');    % default
Kn = trapControlDesign(p, 'designModel', 'double');   % uncompensated

L  = simulateTrapLevitation(p, Kc, x0, tEnd, 'delay', false, ...
                            'Vtia', 0, 'Vmax', 120);
```

`simulateTrapLevitation` options: `'delay'` (acquisition delay on/off),
`'Vtia'` (TIA bias), `'Vmax'` (voltage clamp), `'ref'` (setpoint),
`'trueState'` (bypass the observer and feed the regulator the exact state —
a diagnostic, not implementable).

---

## Verification

```matlab
verifyTrapDeviceGeometry        % geometry vs COMSOL, R^2 >= 0.993, 20 terminals
verifyTrapSensing               % Y-delta reduction vs exact Schur complement
verifyTrapCapToHeight           % capacitance -> height round trip
verifyForceToVoltageRoundTrip   % force -> voltage -> force, machine precision
verifyTrapControl               % controllability, observability, placement
verifyTrapPipeline              % end-to-end checks on every stage
verifyTrapClosedLoop            % 21-state numerical Jacobian of the real loop
```

`verifyTrapControl` reports one expected **FAIL** — the loop delay against the
unstable pole. That is a real, known limitation (see below), not a broken
test.

Use `verifyTrapClosedLoop` rather than `verifyTrapControl` when changing the
actuator layer. `verifyTrapControl` checks the loop against the *design*
model, which is stable by construction and cannot see defects that live in
the voltage stage.

---

## Known limitations

- **Acquisition delay is the binding constraint.** The measured budget is
  roughly 200 µs of total staleness, far tighter than the `p*td < 1` bound
  (874 µs) would suggest — that bound is necessary, not sufficient. With
  `tI2C = 500 µs` and `dtChannel = Tisr/4` the loop does not survive. See
  `presentation/delayBudget.m`.
- **Lateral motion is uncontrolled and leaks into the height estimate.**
  `rx`/`ry` are not observable in the four sensing channels, but they change
  the plate/electrode overlap, so a lateral offset appears as an apparent
  height change of roughly 9 nm per µm. From a ±40 µm start this sets a
  ~1 µm floor on heave accuracy that no bandwidth removes.
- **The lateral channel of the capacitance model is not validated.** Vertical
  and tilt match COMSOL at R² ≥ 0.993, but all Fx/Fy/tauz simulation data in
  this project is corrupted by meshing, so the lateral terms rest on the
  geometric derivation alone.
- **Fringing is deliberately not modelled**, so capacitances are
  systematically low by ~4% RMS. See the header of `trapCapacitance.m`.

---

## Repository conventions

Figures (`*.png`, `*.mp4`) and regenerable caches are gitignored — every plot
script rewrites them. To keep a specific figure, force-add it:

```
git add -f presentation/SomeFigure.png
```

COMSOL exports (`presentation/comsol_*.txt`) and parsed sweep `.mat` files
**are** tracked, because those runs are expensive to reproduce.

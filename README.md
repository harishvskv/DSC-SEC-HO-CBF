# DSC-SEC-HO-CBF

**Decentralised Safety-Critical Secure Edge-Control for Cyber-Physical Microgrids using High-Order Control Barrier Functions.**

A MATLAB implementation of a safety-filtered edge-control architecture for a nonlinear cyber-physical microgrid operating over an unreliable communication link. The repository contains the plant model, the barrier-certificate machinery, the delay-aware event-triggering mechanism, the controller bank, a verification suite, and the numerical studies that exercise all of it.

No proprietary toolbox is required: there is no YALMIP, no MOSEK and no CVX anywhere in the code. The quadratic programs are solved with `quadprog` (Optimization Toolbox), with a dense active-set fallback for installations without it.

---

## 1. Scope of the model

The plant is **nonlinear**. It is never linearised anywhere in this repository, and no small-signal surrogate is substituted for it in any study. The nonlinearities retained in the vector field are

- governor dead-band (GDB) on the speed-governor input,
- generation rate constraint (GRC) on the turbine power derivative,
- composite frequency- and voltage-dependent load with non-unity exponent `gamma`,
- state-of-charge dependent storage derating `Phi(S)`,
- actuator saturation and rate limits on both edge channels.

Each area carries the nine-state vector

```
x = [ df  dPtie  dPm  dPv  z  dPB  S  dPA  dTin ]'   in R^9
```

Under-frequency load shedding is modelled as a staged relay **external** to the vector field (`ufls_relay.m`), so that the drift remains continuous where the barrier analysis requires it.

Two plant configurations are supported: a single isolated area, and a three-area interconnected system with tie-line coupling.

---

## 2. Requirements

| Item | Notes |
|---|---|
| MATLAB | R2023b or later |
| Optimization Toolbox | required, for `quadprog` |
| Statistics and Machine Learning Toolbox | optional, exact binomial confidence intervals |

Developed and run on Windows 11, Intel i7-12700, 15.6 GB RAM. Everything in `src/` is plain `.m` code.

---

## 3. Quick start

```matlab
addpath(genpath('src'));
S = run_all();                      % verification gates, then all studies
make_figures(S);                    % writes figures/ (one plot per file)
```

Partial execution:

```matlab
S = run_all('verify');                              % Layers 0-6 only
S = run_all('studies', struct('S', S, 'n_mc', 50)); % Layer 7, reusing S
S = run_all('all', struct('useCache', true));       % reuse cached constants
```

`run_all` halts at the **first** failed gate and returns `S.gate_failed` naming it. Each layer consumes the outputs of the previous one, so running a study against a stale Lipschitz or contraction-rate struct produces numbers that look plausible and are wrong. The driver makes that dependency order explicit; do not call the study scripts out of order.

Indicative wall-clock times are printed in the header of `run_all.m`. A complete cold run, including the Lipschitz sweep, is of the order of fifty minutes.

---

## 4. Repository layout

```
src/        all MATLAB sources (39 files)
figures/    populated by make_figures.m
results/    populated by run_all.m
```

### Layer 0 — specification and parameters
| File | Purpose |
|---|---|
| `cpmg_spec.m` | single source of truth: state indices, safety limits, unit conversions, under-frequency relay stages |
| `params_isolated.m` | single-area isolated microgrid parameters |
| `params_39bus.m` | three-area parameters and tie-line matrix |

All parameter structs carry a `params_version` stamp. Every simulation entry point checks it and refuses to run against a mismatched struct, which prevents a stale workspace from silently contaminating a study.

### Layer 1 — nonlinear plant
`cpmg_dynamics.m`, `cpmg_integrate.m`, `cpmg_dynamics_multiarea.m`, `cpmg_integrate_multiarea.m`, `ufls_relay.m`

The multi-area drift is assembled block-wise and the input matrix is block diagonal, which is what makes the control architecture genuinely decentralised: one area's inputs cannot reach another's states within an integration step except through the tie-line state.

### Layer 2 — barriers and allowable delay
`barrier_functions.m`, `compute_hcrit.m`, `compute_scatd.m`, `time_to_floor.m`, `validate_scatd.m`

Relative-degree-two barriers, the feasibility floor `h_min`, the critical initial margin `h_crit`, and the safety-critical allowable time delay computed from a finite-time quadratic envelope.

### Layer 3 — cyber layer
`dos_generator.m`, `delay_buffer.m`, `da_etm.m`

`da_etm.m` implements the delay-aware event-triggering mechanism with its three transmission paths (predictive, packet-loss, age-of-data) and the instrumented reconnection logic. Age of data is tracked separately from attack duration; the two are not interchangeable.

### Layer 4 — controllers
`ctrl_hocbf_qp.m`, `ctrl_pi_agc.m`, `ctrl_bank.m`, `estimate_disturbance.m`

`ctrl_bank.m` exposes seven arms for controlled comparison:

| Arm | Description |
|---|---|
| A0 | centralised PI-AGC, no edge layer |
| A1 | PI-AGC with edge fallback, no barrier |
| A2 | HO-CBF filter, no delay awareness |
| A3 | HO-CBF with delay-aware triggering |
| A4 | full proposed scheme |
| A5 | periodic-trigger comparator |
| A6 | delay-compensated predictive comparator |

A6 compensates the **central** path only: the predicted state drives the AGC command, while the edge command is computed from the locally measured state. Driving the edge law from a predicted state produces a comparator that is worse than doing nothing, which is not a meaningful baseline.

Transfer between the central and edge reference laws is bumpless; the shared plant integrator is never reset by either law.

### Layer 5 — simulation engines and metrics
`simulate_case.m`, `simulate_case_multiarea.m`, `metrics.m`, `metrics_multiarea.m`

Three outcomes are recorded separately and never conflated: constraint **violation**, **guarantee-void** (the certified delay bound was exceeded), and **silent failure** (the certificate held but the margin collapsed).

### Layer 6 — verification gates
`verify_plant.m` (29 tests), `verify_lie.m` (24), `verify_etm.m` (12), `verify_qp.m` (14) — 79 assertions in total, covering dimensional consistency, Jacobian-versus-secant Lipschitz bounds, relative-degree claims, trigger exclusivity, packet visibility, QP feasibility across the four-stage cascade, and storage recovery.

The gates are not decoration. Several closed-loop defects in this code base were found only because a gate was added after a trace exposed them, and each threshold involved was individually reasonable and jointly unreachable in the closed loop's actual operating states.

### Layer 7 — studies
`run_isolated.m`, `run_ablation.m`, `run_montecarlo.m`, `run_feasibility_study.m`, `run_39bus.m`, `stats_report.m`

### Layer 8 — constant estimation
`estimate_lipschitz.m`, `estimate_contraction_rates.m`

### Layer 9 — driver, figures, diagnostics
`run_all.m`, `make_figures.m`, `diagnose_case.m`

`make_figures.m` writes one plot per file, no subplots, as 600 dpi TIFF plus a vector PDF, sized to a 3.50 in single-column width and set in Times New Roman.

---

## 5. Reproducibility

- All stochastic studies are seeded. `run_montecarlo.m` sets the generator explicitly at entry, so a repeated run reproduces the same trial set.
- The integrator is fixed-step RK4 with the step declared in `cpmg_spec.m`. No variable-step solver is used anywhere, so results do not depend on tolerance settings.
- Initial conditions, disturbance profiles and attack windows are defined in the study scripts, not hard-coded inside the engines.
- Result archives are **not** committed. `run_all` regenerates them into `results/`.

---

## 6. Contact

V.S.K.V. Harish, Department of Electrical Engineering, Netaji Subhas University of Technology, Dwarka, New Delhi, India.

## 7. License

BSD 3-Clause. Copyright (c) 2026, V.S.K.V. Harish. See `LICENSE`.

Use, modification and redistribution are permitted, including commercially, provided the copyright notice and disclaimer are retained. The third clause additionally bars use of the copyright holder's name to endorse derived work without prior written permission.

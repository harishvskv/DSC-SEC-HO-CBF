function p = params_isolated(spec)
%PARAMS_ISOLATED  Physical, thermal, storage and cyber parameters for the
%   isolated nonlinear cyber-physical microgrid (Case Study 1).
%
%   p = PARAMS_ISOLATED(spec) returns a parameter struct consumed by
%   cpmg_dynamics.m.  Every field printed in Appendix A of the manuscript is
%   taken verbatim from this file, so the paper and the code cannot diverge.
%
%   UNIT DISCIPLINE
%   ---------------
%   Two unit errors present in the legacy implementation are corrected here
%   and are flagged in the comments:
%
%   (E1) BESS energy.  SoC dynamics require ENERGY IN SECONDS OF RATED POWER,
%        not kWh.  Using dSB_dot = -dPB/C_bat_kWh understates the SoC swing by
%        a factor of 3600*S_base/C_bat.  Corrected via E_B_s below.
%
%   (E2) Building thermal capacitance.  R_th*C_th with C_th in kWh/degC gives
%        a time constant in HOURS, but the simulation integrates in SECONDS.
%        The legacy code produced a 16 SECOND building time constant instead
%        of 16 hours, i.e. wrong by 3600x.  Corrected by specifying the
%        aggregate thermal time constant directly in hours.
%
%   See also CPMG_SPEC, CPMG_DYNAMICS, PARAMS_39BUS.

if nargin < 1 || isempty(spec), spec = cpmg_spec(); end
p = struct();
p.name   = 'Isolated nonlinear CPMG';

% Version stamp.  Bumped whenever a FIELD is added or renamed, so that a
% struct left in the workspace from an earlier session is detected at the
% entry point instead of failing deep inside a controller call.
p.params_version = 'params-2.6';
p.spec   = spec;
p.nAreas = 1;
p.isolated = true;

%% ================================================================== %%
%  1. Base quantities
%% ================================================================== %%
p.S_base_MW  = 10.0;                 % microgrid power base            [MW]
p.S_base_kW  = p.S_base_MW*1e3;
p.f_base_Hz  = spec.f_base_Hz;

%% ================================================================== %%
%  2. Rotational and load dynamics
%% ================================================================== %%
p.M      = 8.0;      % equivalent inertia constant                [p.u.-s]
p.D      = 1.2;      % load damping coefficient                   [p.u./p.u.]
p.gamma  = 1.5;      % composite-load frequency-dependence exponent  [-]
p.PL0    = 0.90;     % nominal active load                        [p.u.]
p.Pref   = 0.85;     % scheduled generation set point             [p.u.]

% Composite load model:  P_L(df) = (PL0 + dPL) * (1 + D*df)^gamma
% gamma > 1 makes the load recovery nonlinear in df.  gamma = 0 recovers the
% conventional linear D*df damping term.

%% ================================================================== %%
%  3. Micro-turbine: governor, deadband, turbine, GRC
%% ================================================================== %%
p.Tg     = 0.20;     % governor time constant                     [s]
p.Tt     = 0.30;     % turbine time constant                      [s]
p.R      = 0.04;     % speed droop                                [p.u./p.u.]

% Governor dead-band, Eq.(2).  IEGC / IEEE practice for thermal units is a
% dead-band of the order of +/- 0.03 Hz on a 50 Hz base.
p.df_db      = 0.03/p.f_base_Hz;   % = 6.0e-4 p.u.  (+/- 0.03 Hz)
p.db_smooth_k = 2e3;               % steepness of the smooth surrogate

% Generation rate constraint, Eq.(5)-(6).  Expressed as a RATE limit on the
% mechanical power, not a magnitude limit.  0.1 p.u./min is the conventional
% value for reheat thermal units; micro-turbines are faster.
p.GRC_pu_per_min = 20.0;
p.GRC = p.GRC_pu_per_min/60;       % [p.u./s]
p.grc_smooth_k = 1.0;              % smoothing gain for the tanh surrogate

% Governor mechanical limits (magnitude)
p.Pv_min = -0.50;   p.Pv_max = 0.50;    % valve travel        [p.u.]
p.aw_band = 0.02;   % anti-windup band width [p.u.]; see clampRate.
                    % Set small enough that the physical limit is
                    % tight, large enough that L_F stays bounded.

p.Pm_min = -0.50;   p.Pm_max = 0.50;    % mechanical power    [p.u.]

%% ================================================================== %%
%  4. BESS
%% ================================================================== %%
p.TB       = 0.10;     % BESS power-loop time constant            [s]
p.C_bat_kWh = 500.0;   % installed energy capacity                [kWh]
p.eta_B    = 0.95;     % round-trip efficiency (single-sided)     [-]

% (E1) CORRECTED: energy expressed in seconds of rated power.
%      E_B_s = C_bat[kWh] * 3600 [s/h] / S_base[kW]
p.E_B_s = p.C_bat_kWh*3600/p.S_base_kW;      % = 180 s

% SoC-dependent saturation Phi(dSB), D4.  Symmetric two-sided form:
%   Phi(S) = 0.25*(1+tanh(ks*(S-S_lo)))*(1-tanh(ks*(S-S_hi)))
% Phi -> 1 in the mid-SoC region and Phi -> 0 at BOTH limits, so the BESS
% column of G(x) degenerates at either extreme.  This is the mechanism that
% can render the QP infeasible during prolonged DoS.
p.SoC_nom   = spec.safety.SoC_nom;
p.SoC_min   = spec.safety.SoC_min;
p.SoC_max   = spec.safety.SoC_max;
p.SoC_marg  = 0.05;                                % knee margin
p.SoC_lo    = p.SoC_min + p.SoC_marg;
p.SoC_hi    = p.SoC_max - p.SoC_marg;
p.ks        = 60.0;                                % Phi steepness
p.Phi_min   = 1e-3;                                % degeneracy detection floor

% BESS actuator limits (also imposed as QP bounds)

% Enforced (effective) storage bounds, per spec.safety.SoC_anchor.  The
% barrier is anchored here rather than at the physical bound, because Phi has
% already collapsed by the time the physical bound is reached.
p.SoC_enf_lo = p.SoC_lo;
p.SoC_enf_hi = p.SoC_hi;
phiAt = @(S) 0.25*(1+tanh(p.ks*(S-p.SoC_lo)))*(1-tanh(p.ks*(S-p.SoC_hi)));
p.Phi_at_anchor = min(phiAt(p.SoC_enf_lo), phiAt(p.SoC_enf_hi));
assert(p.Phi_at_anchor >= spec.safety.Phi_floor_at_anchor, ...
    ['params: Phi at the enforced storage bound is %.3f, below the required ' ...
     'floor %.3f. Widen SoC_marg or reduce ks.'], ...
     p.Phi_at_anchor, spec.safety.Phi_floor_at_anchor);

p.uB_min = -0.60;  p.uB_max = 0.60;                % [p.u.]
p.uB_rate  = 12.0;    % [p.u./s] BESS command slew limit, QP rate row.
% Corresponds to a full-range traverse (uB_min to uB_max, 1.2 p.u.) in
% 100 ms, which is within the capability of a grid-scale four-quadrant
% inverter and consistent with fast frequency response requirements.
% The earlier value 3.0 p.u./s needed 200 ms for full output, longer than
% the time in which the barrier is lost from a thin margin, so the
% actuator could not physically arrest the fall.
p.uB_full_traverse_s = (p.uB_max - p.uB_min)/p.uB_rate;

%% ================================================================== %%
%  5. AIAC fleet and building thermal model
%% ================================================================== %%
p.Tac    = 0.20;     % aggregate compressor / inverter time constant  [s]
p.COP    = 2.50;     % nominal coefficient of performance             [-]

% (E2) CORRECTED: aggregate thermal dynamics specified directly, in a form
%      that is dimensionally safe when time is integrated in SECONDS.
%
%      dTin_dot = a_T*(dTout - dTin) - b_T*dPA
%
%      a_T = 1/tau_th   with tau_th in SECONDS
%      b_T = k_th*a_T   so that the steady-state gain is dTin = -k_th*dPA
p.tau_th_hr = 4.0;                        % aggregate building time constant [h]
p.a_T = 1/(p.tau_th_hr*3600);             % [1/s]

% Ambient and comfort set points (needed before k_th can be derived)
p.Tout_nom = 32.0;                        % nominal outdoor temperature [degC]
p.Tin_ref  = spec.safety.Tin_ref_C;       % comfort set point           [degC]
p.PA_nom   = 0.20;                        % nominal AIAC fleet consumption [p.u.]

% Thermal gain is DERIVED, not assumed, so that the nominal operating point
% is an exact equilibrium of the thermal state:
%     0 = a_T*(Tout_nom - Tin_ref) - b_T*PA_nom
% Choosing k_th independently (as the legacy model did) leaves a constant
% residual and the indoor temperature drifts even with no control action.
p.k_th_pu = (p.Tout_nom - p.Tin_ref)/p.PA_nom;   % [degC per p.u.]
p.b_T     = p.k_th_pu*p.a_T;                     % [degC/(s*p.u.)]

% SI-equivalent parameters, reported in Appendix A for interpretability
p.Rth_eq_CperkW  = p.k_th_pu/(p.COP*p.S_base_kW);       % [degC/kW]
p.Cth_eq_kWhperC = p.tau_th_hr/p.Rth_eq_CperkW;         % [kWh/degC]

% CONSEQUENCE OF THE (E2) CORRECTION, stated here so it is not discovered
% later: with a physically correct 4 h aggregate thermal time constant, the
% indoor temperature moves by O(1e-3) degC over a 3 s DoS event.  The thermal
% comfort barrier is therefore INACTIVE for short attacks and becomes active
% only for sustained (minutes-scale) events.  Any claim of tight transient
% comfort regulation over a 3 s window is an artefact of the legacy 3600x
% unit error and must not be repeated in the manuscript.

% Optional COP derating with temperature difference (spec.nl.cop_temperature)
p.COP_dTnom = 10.0;      % reference (Tout - Tin) at which COP = p.COP  [degC]
p.COP_slope = 0.02;      % fractional COP loss per degC above reference [1/degC]
p.COP_min   = 1.5;  p.COP_max = 4.0;

% Comfort band
p.Tin_min   = spec.safety.Tin_min_C;
p.Tin_max   = spec.safety.Tin_max_C;

% AIAC actuator limits.  Hui2019 reports an explicit compressor rate limiter;
% it is imposed as a QP RATE ROW rather than inside the plant, because a rate
% limiter in the plant would destroy control-affineness in u and invalidate
% the HO-CBF construction.
p.uA_min = -0.30;  p.uA_max = 0.30;       % [p.u.] about the AIAC operating point
p.uA_rate  = 0.10;    % [p.u./s] AIAC compressor slew limit, after Hui2019.
% CONSEQUENCE, to be stated rather than discovered: at this rate the AIAC
% moves 0.001 p.u. per 10 ms controller period, so it contributes
% essentially nothing over the sub-second window in which the frequency
% barrier is at risk. Fast frequency safety therefore rests entirely on
% the BESS, which makes the Phi(SoC) degeneracy more consequential, not
% less. The AIAC provides slow support and comfort-bounded energy relief.
p.uA_full_traverse_s = (p.uA_max - p.uA_min)/p.uA_rate;

%% ================================================================== %%
%  6. Area control error
%% ================================================================== %%
% Isolated topology: no tie line, so ACE = beta*df.
p.beta = 1/p.R + p.D;                     % natural frequency response [p.u./p.u.]
p.T_ij = 0;                               % no synchronising coefficients

%% ================================================================== %%
%  7. Central AGC (WAN side, subject to DoS)
%% ================================================================== %%
% Central AGC to the governor: proportional plus integral on ACE.
p.agc.Kp = 0.40;
p.agc.Ki = 0.50;
% RETUNED from 0.10.  With integral action removed from the edge reference
% (Ki_e = 0), restoration falls entirely to the AGC, and Ki = 0.10 left it
% taking 26.6 s to return inside the IEGC band because the edge integral
% had quietly been doing the AGC's job.  A gain sweep on the closed loop
% gives settling times of 26.6, 5.0, 2.7 and 1.3 s at Ki = 0.10, 0.50,
% 0.80 and 1.20, with no overshoot up to 0.80.  Oscillation begins near
% Ki = 1.2 (3 sign changes) and the loop is clearly unstable by Ki = 5.
% Ki = 0.50 therefore sits roughly a factor of 2.4 below the onset of
% oscillation and a factor of 10 below instability.  The nadir is
% insensitive to this gain (-0.5742 to -0.5734 Hz across the sweep):
% it is set by the transient, which the barrier governs, not by
% restoration speed.
p.agc.u_min = -0.50;  p.agc.u_max = 0.50;

% Edge dispatch reference.  NO INTEGRAL TERM: zACE settles at a nonzero
% constant after a load step, so integral action on the edge channels
% commands the storage to supply steady-state energy indefinitely.
% See ctrl_pi_agc for the trace that exposed this.
p.agc.Kp_e  = 0.50;
p.agc.Ki_e  = 0.00;
p.agc.K_soc = 0.50;   % BESS charge restoration gain [p.u. per SoC]
% restoration time constant approx E_B/(eta_B*K_soc)
p.agc.T_soc_restore = p.E_B_s/(p.eta_B*max(p.agc.K_soc,eps));

%% ================================================================== %%
%  8. Cyber layer
%% ================================================================== %%
p.cyber.T_pkt      = 0.050;   % WAN packet period, central AGC to edge [s]
p.cyber.tau_pred_max = 1.000;  % cap on the A6 prediction horizon [s].
                              % Beyond this the frozen-input prediction
                              % is extrapolation, and an uncapped horizon
                              % would flatter the proposed method.
p.cyber.tau_nom_lo = 0.020;   % nominal WAN latency, lower  [s]
p.cyber.tau_nom_hi = 0.040;   % nominal WAN latency, upper  [s]

% WORST HEALTHY AGE OF DATA.  Must sit strictly below the delay
% feasibility limit tau_max from Layer 2 (0.1463 s for this plant).
% The earlier value tau_nom_hi = 0.100 gave a worst healthy age of
% 0.150 s, ABOVE tau_max, so the tightened HO-CBF constraint would have
% been infeasible during normal operation with no attack at all.
p.cyber.age_healthy_max = p.cyber.T_pkt + p.cyber.tau_nom_hi;

% Worst-case ATTACK DURATION simulated anywhere in the study.  The Monte
% Carlo sweep draws durations up to 3.4 s, so this must exceed that.
p.cyber.dos_dur_max = 3.500;

% AGE OF DATA IS NOT ATTACK DURATION.  After a blackout of length D the
% edge is holding a packet generated one period and one flight time
% BEFORE the attack began, and the first post-attack packet does not
% arrive until one more period and flight time AFTER it ends:
%
%     age_dos_max = D + T_pkt + tau_nom_hi
%
% Conflating the two is what sized the buffer wrongly.  It is also why a
% manuscript cannot describe a 3 s blackout as producing a 2 s data age.
p.cyber.age_dos_max = p.cyber.dos_dur_max + p.cyber.T_pkt + p.cyber.tau_nom_hi;
p.cyber.tau_dos_hi  = p.cyber.age_dos_max;   % retained name, same quantity

% Delay-buffer horizon, with 50 percent headroom over the worst age.
p.cyber.tau_buffer = 1.5*p.cyber.age_dos_max;
assert(p.cyber.tau_buffer > p.cyber.age_dos_max, ...
    'params: delay buffer horizon must exceed the worst-case age of data');
assert(p.cyber.age_healthy_max < p.cyber.tau_nom_hi + p.cyber.T_pkt + eps, ...
    'params: healthy age of data is inconsistent with T_pkt and tau_nom_hi');

% Event-trigger settings (Layer 5)
p.etm.tau_trigger = 0.120;    % age-of-data backstop; healthy_max < this < tau_max
p.etm.n_miss      = 2;        % consecutive missed packets that declare a DoS
p.etm.T_react     = 0.150;    % predictive-trigger horizon [s]
p.etm.T_dwell     = 1.000;    % minimum isolation dwell, anti-chatter [s]
p.etm.n_ok        = 3;        % fresh packets required before reconnection
p.etm.h_recover   = 0.70;     % fraction of the safe set required to reconnect
p.etm.k_recover   = 3.0;      % T_pred must exceed k_recover*T_react to
                              % reconnect. Was hard-coded at 3 inside da_etm.
p.etm.T_rocof     = 0.020;    % RoCoF measurement filter lag [s]; adds to the
                              % timing budget alongside T_pkt and WAN latency
p.etm.T_bump      = 2.000;    % bumpless-transfer offset decay [s].
                              % Without it, switching between two reference
                              % laws that share the plant integrator zACE
                              % steps the command and the resulting RoCoF
                              % spike re-fires the predictive trigger.
assert(p.cyber.age_healthy_max < p.etm.tau_trigger, ...
    'params: trigger would fire on healthy jitter');

%% ================================================================== %%
%  9. Initial condition (equilibrium at nominal load, u = 0)
%% ================================================================== %%
ix = spec.ix;
x0 = zeros(spec.nx,1);
x0(ix.df)    = 0.0;
x0(ix.dPtie) = 0.0;
x0(ix.dPm)   = 0.0;
x0(ix.dPv)   = 0.0;
x0(ix.zACE)  = 0.0;
x0(ix.dPB)   = 0.0;
x0(ix.dSB)   = p.SoC_nom;      % ABSOLUTE SoC, not a deviation
x0(ix.dPA)   = 0.0;
x0(ix.dTin)  = 0.0;            % deviation from Tin_ref
p.x0 = x0;

end

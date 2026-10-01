function spec = cpmg_spec()
%CPMG_SPEC  Frozen specification for the DSC-SEC-HOCBF framework.
%
%   spec = CPMG_SPEC() returns the single authoritative definition of the
%   canonical state vector, safety tiers, sign conventions and units used by
%   every other file in this repository.  NO other file may hard-code a state
%   index, a safety threshold or a sign convention.  If a definition needs to
%   change, it changes here and only here.
%
%   This file exists to eliminate the class of defect in which the manuscript,
%   the controller and the plant disagree about what x(4) means or what the
%   safety limit is.
%
%   ---------------------------------------------------------------------
%   D1  CANONICAL STATE VECTOR (identical for isolated and multi-area)
%   ---------------------------------------------------------------------
%       x_i = [ df_i ; dPtie_i ; dPm_i ; dPv_i ; zACE_i ;
%               dPB_i ; dSB_i ; dPA_i ; dTin_i ]  in R^9
%
%   For an ISOLATED microgrid the dimension is NOT reduced.  Instead
%   dPtie_i is held identically zero by setting the tie-line coupling input
%   to zero, and ACE_i degenerates to beta_i*df_i.  This keeps one index map
%   across all topologies.
%
%   ---------------------------------------------------------------------
%   D2  SAFETY TIERS (Indian Electricity Grid Code 2023 / NLDC OP)
%   ---------------------------------------------------------------------
%   Reference frequency 50.000 Hz.  IEGC normal operating band is
%   49.900 - 50.050 Hz.  Default under-frequency relay (UFR) stages are
%   49.40 / 49.20 / 49.00 / 48.80 Hz.
%
%   The CONTROL BARRIER FUNCTION is anchored to UFR Stage-1 (49.40 Hz):
%
%       h_f(x) = df - df_min ,      df_min = -0.012 p.u.  (-0.60 Hz)
%
%   The IEGC normal band (-0.0020 / +0.0010 p.u.) is a PERFORMANCE metric
%   only.  It is never referred to as a safety constraint anywhere in the
%   code or the manuscript.
%
%   ---------------------------------------------------------------------
%   D3  SIGN CONVENTIONS  (after Hui et al., IEEE TIE 66(2) 2019; 68(3) 2021)
%   ---------------------------------------------------------------------
%   dPB   > 0  : BESS DISCHARGING, i.e. injecting power into the grid.
%                Enters the swing equation with a PLUS sign.
%                SoC therefore falls:  dSB_dot = -dPB / E_B.
%
%   dPA   > 0  : AIAC fleet CONSUMING MORE electrical power (more cooling).
%                Enters the swing equation with a MINUS sign.
%                Indoor temperature therefore falls.
%                Consequently LgLf_h for the AIAC channel is NEGATIVE and a
%                frequency-supporting action requires u_A < 0 (load shedding).
%
%   Hui2021 Eq.(14) writes +dP_IAC because their variable is the regulation
%   CAPACITY PROVIDED (the negative of consumption).  Both conventions are
%   self-consistent; mixing them is not.  This repository uses CONSUMPTION
%   throughout, because the ETP thermal coupling is written in terms of the
%   electrical power actually drawn by the compressors.
%
%   dPtie > 0  : NET EXPORT from the area.  Enters the swing with a MINUS.
%   dPm   > 0  : mechanical power above the scheduled set point.
%
%   ---------------------------------------------------------------------
%   D4  BESS SoC-DEPENDENT SATURATION  Phi(dSB)  [CONFIRMED, retained]
%   ---------------------------------------------------------------------
%   Phi multiplies the BESS control channel in G(x) and is therefore present
%   in LgLf_h.  Phi -> 0 at BOTH SoC limits, which degenerates the relative
%   degree and can render the QP infeasible during prolonged attacks.  This
%   is a property of the physics, not a defect, and it is analysed explicitly
%   rather than hidden.
%
%   ---------------------------------------------------------------------
%   D5  CLF RETAINED
%   ---------------------------------------------------------------------
%   The QP carries a genuine CLF row WITH a slack variable, uniformly in the
%   isolated and multi-area cases.  Declared here so no downstream file can
%   quietly omit it.
%
%   Author : DSC-SEC-HOCBF project
%   Repo   : https://github.com/harishvskv/DSC-SEC-HO-CBF

%% ------------------------------------------------------------------ %%
%  Version stamp (bump on any change; logged into every results file)
%% ------------------------------------------------------------------ %%
spec.version      = 'spec-2.0';
spec.description  = 'Canonical R^9 CPMG spec, IEGC-2023 safety tiers, Hui sign convention';

%% ------------------------------------------------------------------ %%
%  D1  State index map
%% ------------------------------------------------------------------ %%
spec.nx = 9;                 % states per area
spec.nu = 2;                 % edge control inputs per area [u_B ; u_A]
spec.nd = 3;                 % disturbance channels [dP_res ; dP_L ; dT_out]

spec.ix.df    = 1;           % frequency deviation                     [p.u.]
spec.ix.dPtie = 2;           % net tie-line export deviation           [p.u.]
spec.ix.dPm   = 3;           % turbine mechanical power deviation      [p.u.]
spec.ix.dPv   = 4;           % governor valve position deviation       [p.u.]
spec.ix.zACE  = 5;           % integral of ACE                         [p.u.-s]
spec.ix.dPB   = 6;           % BESS power injection deviation          [p.u.]
spec.ix.dSB   = 7;           % BESS state of charge, ABSOLUTE          [fraction]
spec.ix.dPA   = 8;           % AIAC aggregate CONSUMPTION deviation    [p.u.]
spec.ix.dTin  = 9;           % indoor temperature deviation from Tin_ref [degC]

% NOTE on x(7): unlike every other entry this is an ABSOLUTE quantity in
% [0,1], not a deviation.  The SoC barrier and the saturation Phi are both
% defined on absolute SoC, so storing a deviation would require carrying the
% operating point separately.  The symbol is kept as dSB for continuity with
% the manuscript notation.
%
% NOTE on x(9) and d(3): both are deviations, x(9) from spec.safety.Tin_ref_C
% and d(3) from p.Tout_nom.  Absolute temperatures are recovered as
% Tin  = Tin_ref  + x(9)  and  Tout = Tout_nom + d(3).

spec.iu.uB    = 1;           % BESS power command                      [p.u.]
spec.iu.uA    = 2;           % AIAC power command                      [p.u.]

spec.id.dPres = 1;           % renewable injection deviation           [p.u.]
spec.id.dPL   = 2;           % load demand deviation                   [p.u.]
spec.id.dTout = 3;           % ambient temperature deviation           [degC]

spec.stateNames = { 'df','dPtie','dPm','dPv','zACE','dPB','dSB','dPA','dTin' };
spec.stateUnits = { 'pu','pu','pu','pu','pu*s','pu','frac','pu','degC' };
spec.inputNames = { 'u_BESS','u_AIAC' };
spec.distNames  = { 'dP_res','dP_L','dT_out' };

% Consistency guard: index map must be a permutation of 1..nx
assert(isequal(sort(cell2mat(struct2cell(spec.ix))'),1:spec.nx), ...
    'cpmg_spec: state index map is not a permutation of 1..nx');

%% ------------------------------------------------------------------ %%
%  D2  Frequency tiers, IEGC 2023
%% ------------------------------------------------------------------ %%
spec.f_base_Hz = 50.0;                       % system base frequency

hz2pu = @(fHz) (fHz - spec.f_base_Hz)/spec.f_base_Hz;
spec.hz2pu = hz2pu;
spec.pu2hz = @(pu) pu*spec.f_base_Hz;        % deviation in pu -> deviation in Hz

% IEGC 2023 normal operating band (PERFORMANCE metric, not a safety set)
spec.band.f_lo_Hz   = 49.900;
spec.band.f_hi_Hz   = 50.050;
spec.band.df_lo_pu  = hz2pu(49.900);         % -0.0020 p.u.
spec.band.df_hi_pu  = hz2pu(50.050);         % +0.0010 p.u.
spec.band.label     = 'IEGC-2023 normal operating band (performance only)';

% Under-frequency relay stages, NLDC default settings
spec.ufls.f_Hz      = [49.40 , 49.20 , 49.00 , 48.80];
spec.ufls.df_pu     = hz2pu(spec.ufls.f_Hz); % [-0.0120 -0.0160 -0.0200 -0.0240]
spec.ufls.shed_frac = [0.05  , 0.05  , 0.05  , 0.05 ];  % load shed per stage
spec.ufls.delay_s   = [0.20  , 0.20  , 0.20  , 0.20 ];  % relay operating time
spec.ufls.reset_Hz  = 49.90;                 % hysteresis reset level
spec.ufls.stageName = {'Stage-1','Stage-2','Stage-3','Stage-4'};

% >>> THE safety threshold used by the control barrier function <<<
spec.safety.stage_index = 1;                 % UFR Stage-1
spec.safety.df_min_pu   = spec.ufls.df_pu(spec.safety.stage_index);   % -0.0120
spec.safety.df_min_Hz   = spec.ufls.f_Hz(spec.safety.stage_index);    %  49.40
spec.safety.source      = 'IEGC 2023 / NLDC Operating Procedure, UFR Stage-1 (49.40 Hz)';

% Over-frequency limit.  Not standard-anchored; declared as an engineering
% choice and reported as such.  Only used for reporting, not for the CBF.
spec.safety.df_max_pu   = hz2pu(50.50);
spec.safety.df_max_note = 'Engineering choice (50.50 Hz). Not IEGC-mandated. Reporting only.';

%% ------------------------------------------------------------------ %%
%  Auxiliary safe sets enforced as CBF rows in the QP
%% ------------------------------------------------------------------ %%
spec.safety.SoC_min   = 0.20;                % BESS PHYSICAL lower SoC bound
spec.safety.SoC_max   = 0.90;                % BESS PHYSICAL upper SoC bound
spec.safety.SoC_nom   = 0.50;                % nominal / reference SoC

% ENFORCED versus PHYSICAL storage bounds
% ---------------------------------------
% The BESS control channel is scaled by Phi(SoC), which collapses over a
% margin band adjacent to each physical bound.  Anchoring the storage barrier
% at the PHYSICAL limit would ask the barrier to defend a region in which the
% actuator has no authority: between SoC_min and SoC_min + margin the barrier
% reports a positive margin while Phi is already near zero.
%
% The barrier is therefore anchored at the EFFECTIVE limits SoC_lo and SoC_hi
% defined in the parameter files, where Phi is still one half of its nominal
% value.  This is the largest storage set the actuator can actually render
% forward invariant.  The physical bounds remain in force as reported limits
% and as the domain of Phi; they are simply not the set the QP enforces.
%
% Consequence to state in the manuscript: the enforceable safe set is strictly
% smaller than the physical one, and the gap is set by the derating margin and
% by the Phi steepness.  This is the constructive form of Proposition 1.
spec.safety.SoC_anchor = 'effective';        % 'effective' | 'physical'
spec.safety.Phi_floor_at_anchor = 0.45;      % required Phi at the enforced bound

spec.safety.Tin_min_C = 20.0;                % ASHRAE thermal comfort lower bound
spec.safety.Tin_max_C = 24.0;                % ASHRAE thermal comfort upper bound
spec.safety.Tin_ref_C = 22.0;                % comfort set point

%% ------------------------------------------------------------------ %%
%  D3  Sign conventions, recorded so downstream code can assert on them
%% ------------------------------------------------------------------ %%
spec.sign.swing_dPB   = +1;   % BESS injection adds to power balance
spec.sign.swing_dPA   = -1;   % AIAC consumption subtracts from power balance
spec.sign.swing_dPtie = -1;   % net export subtracts from power balance
spec.sign.thermal_dPA = -1;   % more AIAC power lowers indoor temperature
spec.sign.soc_dPB     = -1;   % discharging lowers SoC
spec.sign.reference   = 'Hui et al., IEEE TIE 66(2):1413, 2019; IEEE TIE 68(3):2725, 2021';

%% ------------------------------------------------------------------ %%
%  D4 / D5  Structural switches
%% ------------------------------------------------------------------ %%
spec.use_Phi_in_G = true;     % D4: Phi(dSB) present in G(x) and hence in LgLf_h
spec.use_CLF      = true;     % D5: CLF row with slack retained in the QP

%% ------------------------------------------------------------------ %%
%  Nonlinearity switches for the plant (see cpmg_dynamics.m)
%% ------------------------------------------------------------------ %%
spec.nl.governor_deadband = true;    % Section II-B
spec.nl.turbine_GRC       = true;    % Section II-B
spec.nl.load_frequency    = true;    % gamma-exponent composite load
spec.nl.bess_Phi          = true;    % Section II-C
spec.nl.cop_temperature   = false;   % optional COP(T_out,T_in) derating

% Deadband representation used INSIDE THE PLANT.  'hard' is the physical
% dead-zone of Eq.(2).  The barrier module never differentiates through it:
% the deadband acts on dPv, which first appears at the THIRD Lie derivative
% of h_f, whereas the HO-CBF has relative degree 2.  A smooth surrogate is
% therefore unnecessary for the controller and is provided only for analysis.
spec.nl.deadband_mode = 'hard';      % 'hard' | 'smooth'
spec.nl.grc_mode      = 'hard';      % 'hard' | 'smooth'

%% ------------------------------------------------------------------ %%
%  Numerical settings
%% ------------------------------------------------------------------ %%
spec.num.Ts_plant   = 1e-3;   % plant integration step [s] (RK4)
spec.num.Ts_ctrl    = 1e-2;   % edge controller / ZOH sampling period [s]
spec.num.integrator = 'rk4';  % 'rk4' | 'euler'
spec.num.eps_guard  = 1e-9;   % guard against division / fractional powers of <=0

%% ------------------------------------------------------------------ %%
%  Debug hooks, off by default
%% ------------------------------------------------------------------ %%
%  cpmg_dynamics_multiarea can assert that the stacked input matrix is block
%  diagonal on every call, which is the structural condition that makes the
%  controller decentralised.  It is off by default because the check runs
%  inside the integration loop.  verify_plant T16 asserts the same property
%  once per run, which is the normal route; this hook exists for the case
%  where a violation is suspected mid-trajectory.
spec.debug.checkBlockDiag = false;

assert(mod(round(spec.num.Ts_ctrl/spec.num.Ts_plant),1)==0, ...
    'cpmg_spec: Ts_ctrl must be an integer multiple of Ts_plant');

end

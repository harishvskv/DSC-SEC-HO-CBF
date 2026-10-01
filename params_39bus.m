function P = params_39bus(spec)
%PARAMS_39BUS  Parameters for the IEEE 39-bus, 10-generator, 3-area
%   nonlinear cyber-physical microgrid (Case Study 2).
%
%   P = PARAMS_39BUS(spec) returns a 1x3 struct array of per-area parameter
%   structs, each of the SAME form returned by params_isolated, plus a global
%   field P(1).net holding the network coupling data.
%
%   Every area carries the full set of nonlinearities: governor dead-band,
%   generation rate constraint, gamma-dependent composite load, and the
%   SoC-dependent BESS saturation Phi.  There is no linearised variant of the
%   plant anywhere in this repository.
%
%   TIE-LINE REPRESENTATION
%   -----------------------
%   dPtie_i_dot = 2*pi*sum_j T_ij*(df_i - df_j)
%
%   This is the conventional synchronising-power representation used
%   throughout the LFC literature.  It is a LINEAR coupling and is declared
%   as such in Assumption 3 of the manuscript.  It is the only linear element
%   of the model; the machine, load, storage and thermal dynamics are all
%   nonlinear.  A sin(delta_i - delta_j) coupling would require augmenting the
%   state with rotor angles and is left as stated future work.
%
%   See also CPMG_SPEC, PARAMS_ISOLATED, CPMG_DYNAMICS_MULTIAREA.

if nargin < 1 || isempty(spec), spec = cpmg_spec(); end

nA = 3;

% Per-area heterogeneous data.  Area 3 is deliberately the low-inertia,
% heavily loaded area and is the stress case throughout the paper.
M_a     = [10.0 ,  8.0 ,  6.0 ];      % inertia constant            [p.u.-s]
D_a     = [ 1.0 ,  1.2 ,  0.9 ];      % load damping                [p.u./p.u.]
Tt_a    = [ 0.30,  0.33,  0.35];      % turbine time constant       [s]
Tg_a    = [ 0.08,  0.07,  0.07];      % governor time constant      [s]
R_a     = [ 0.05,  0.05,  0.05];      % droop                       [p.u./p.u.]
gam_a   = [ 1.5 ,  1.5 ,  1.5 ];      % load exponent                [-]
PL0_a   = [ 0.90,  0.72,  1.08];      % nominal load                [p.u.]
Sbase_a = [ 100 ,  100 ,  100 ];      % area power base             [MW]
GRCm_a  = [ 10.0,  10.0,  10.0];      % GRC                    [p.u./min]
Cbat_a  = [5000 , 5000 , 5000 ];      % BESS energy capacity        [kWh]

for i = nA:-1:1
    p = struct();
    p.name       = sprintf('IEEE 39-bus Area %d', i);
    p.params_version = 'params-2.6';
    p.spec       = spec;
    p.areaIndex  = i;
    p.nAreas     = nA;
    p.isolated   = false;

    % ---- base ----
    p.S_base_MW = Sbase_a(i);
    p.S_base_kW = p.S_base_MW*1e3;
    p.f_base_Hz = spec.f_base_Hz;

    % ---- rotational / load ----
    p.M     = M_a(i);
    p.D     = D_a(i);
    p.gamma = gam_a(i);
    p.PL0   = PL0_a(i);
    p.Pref  = PL0_a(i);

    % ---- turbine chain ----
    p.Tg = Tg_a(i);
    p.Tt = Tt_a(i);
    p.R  = R_a(i);

    p.df_db        = 0.03/p.f_base_Hz;      % +/- 0.03 Hz dead-band
    p.db_smooth_k  = 2e3;
    p.GRC_pu_per_min = GRCm_a(i);
    p.GRC          = p.GRC_pu_per_min/60;
    p.grc_smooth_k = 1.0;
    p.Pv_min = -0.50;  p.Pv_max = 0.50;
    p.aw_band = 0.02;   % anti-windup band width [p.u.]; see clampRate.
                        % Set small enough that the physical limit is
                        % tight, large enough that L_F stays bounded.

    p.Pm_min = -0.50;  p.Pm_max = 0.50;

    % ---- BESS ----
    p.TB        = 0.10;
    p.C_bat_kWh = Cbat_a(i);
    p.eta_B     = 0.95;
    p.E_B_s     = p.C_bat_kWh*3600/p.S_base_kW;    % [s] -- see (E1) note

    p.SoC_nom  = spec.safety.SoC_nom;
    p.SoC_min  = spec.safety.SoC_min;
    p.SoC_max  = spec.safety.SoC_max;
    p.SoC_marg = 0.05;
    p.SoC_lo   = p.SoC_min + p.SoC_marg;
    p.SoC_hi   = p.SoC_max - p.SoC_marg;
    p.ks       = 60.0;
    p.Phi_min  = 1e-3;
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
    p.uB_min   = -0.60;  p.uB_max = 0.60;
    p.uB_rate  = 12.0;    % [p.u./s] BESS command slew limit, QP rate row.
    % Corresponds to a full-range traverse (uB_min to uB_max, 1.2 p.u.) in
    % 100 ms, which is within the capability of a grid-scale four-quadrant
    % inverter and consistent with fast frequency response requirements.
    % The earlier value 3.0 p.u./s needed 200 ms for full output, longer than
    % the time in which the barrier is lost from a thin margin, so the
    % actuator could not physically arrest the fall.
    p.uB_full_traverse_s = (p.uB_max - p.uB_min)/p.uB_rate;

    % ---- AIAC / thermal ----
    p.Tac = 0.20;
    p.COP = 2.50;
    p.tau_th_hr = 4.0;
    p.a_T       = 1/(p.tau_th_hr*3600);
    p.Tout_nom  = 32.0;
    p.Tin_ref   = spec.safety.Tin_ref_C;
    p.PA_nom    = 0.20;
    % Derived so that the nominal point is an exact thermal equilibrium
    p.k_th_pu = (p.Tout_nom - p.Tin_ref)/p.PA_nom;
    p.b_T     = p.k_th_pu*p.a_T;
    p.Rth_eq_CperkW  = p.k_th_pu/(p.COP*p.S_base_kW);
    p.Cth_eq_kWhperC = p.tau_th_hr/p.Rth_eq_CperkW;

    p.COP_dTnom = 10.0;  p.COP_slope = 0.02;
    p.COP_min   = 1.5;   p.COP_max   = 4.0;

    p.Tin_min  = spec.safety.Tin_min_C;
    p.Tin_max  = spec.safety.Tin_max_C;
    p.uA_min   = -0.30;  p.uA_max = 0.30;
    p.uA_rate  = 0.10;    % [p.u./s] AIAC compressor slew limit, after Hui2019.
    % CONSEQUENCE, to be stated rather than discovered: at this rate the AIAC
    % moves 0.001 p.u. per 10 ms controller period, so it contributes
    % essentially nothing over the sub-second window in which the frequency
    % barrier is at risk. Fast frequency safety therefore rests entirely on
    % the BESS, which makes the Phi(SoC) degeneracy more consequential, not
    % less. The AIAC provides slow support and comfort-bounded energy relief.
    p.uA_full_traverse_s = (p.uA_max - p.uA_min)/p.uA_rate;

    % ---- ACE ----
    p.beta = 1/p.R + p.D;          % natural frequency bias [p.u./p.u.]

    % ---- central AGC ----
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

    % ---- cyber ----
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
    p.etm.T_rocof     = 0.020;    % RoCoF measurement filter lag [s]
    p.etm.T_bump      = 2.000;    % bumpless-transfer offset decay [s].
                                  % Without it, switching between two reference
                                  % laws that share the plant integrator zACE
                                  % steps the command and the resulting RoCoF
                                  % spike re-fires the predictive trigger.
    assert(p.cyber.age_healthy_max < p.etm.tau_trigger, ...
        'params: trigger would fire on healthy jitter');

    % ---- initial condition ----
    x0 = zeros(spec.nx,1);
    x0(spec.ix.dSB) = p.SoC_nom;
    p.x0 = x0;

    P(i) = p; %#ok<AGROW>
end

%% ------------------------------------------------------------------ %%
%  Network coupling (shared by all areas)
%% ------------------------------------------------------------------ %%
net = struct();
net.nAreas = nA;

% Synchronising coefficients T_ij [p.u./rad].  Symmetric, zero diagonal.
T = [ 0.00 , 0.20 , 0.25 ;
      0.20 , 0.00 , 0.15 ;
      0.25 , 0.15 , 0.00 ];
assert(isequal(T,T'), 'params_39bus: T_ij must be symmetric');
assert(all(diag(T)==0), 'params_39bus: T_ii must be zero');
net.T_ij = T;

% Disturbance distribution factors across areas (geographic dispersion)
net.dist_share = [1.00 , 0.80 , 1.20];

% Global index of a per-area state within the stacked 27-state vector
net.gidx = @(area,local) (area-1)*spec.nx + local;

P(1).net = net;
for i = 2:nA, P(i).net = net; end

end

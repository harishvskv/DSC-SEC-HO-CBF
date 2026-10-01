function [ug, u_edge, cs, ci] = ctrl_bank(arm, t, x, x_wan, rocof, net, p, ctx, cs)
%CTRL_BANK  Dispatcher for the six ablation arms.
%
%   [ug, u_edge, cs, ci] = CTRL_BANK(arm, t, x, x_wan, rocof, net, p, ctx, cs)
%
%   WHY AN ABLATION BANK EXISTS AT ALL
%   ----------------------------------
%   The submitted manuscript compared the proposed framework against a single
%   PI-AGC baseline that had NO battery and NO air conditioning load.  In the
%   legacy code the baseline edge command was literally zeros(6,1).  None of
%   the reported improvement could therefore be attributed to the controller
%   rather than to the two extra actuators.  Reviewe[ug, ue] = ctrl_pi_agc(xp, p, 'agc_edge');r 1 comment 11, Reviewer 2
%   comment 7 and Reviewer 3 comment 3 all say versions of this, and Reviewer 3
%   states the requirement precisely: show whether the gain comes from the
%   switching, the barrier, the tightening, or simply the added authority.
%
%   ARMS
%     'A0'  PI-AGC, no edge actuation.               The legacy baseline.
%     'A1'  PI-AGC with IDENTICAL BESS and AIAC authority.
%           THIS IS THE ARM THAT DECIDES WHETHER THE PAPER HAS A
%           CONTRIBUTION.  If A4 does not beat A1, the improvement was the
%           actuators, not the controller, and that must be reported.
%     'A2'  Safety filter with NO staleness detection: the filter runs on the
%           local state, but the reference it tracks is the last AGC command
%           received, however old.  Isolates the value of the DA-ETM.
%     'A3'  HO-CBF QP with DA-ETM but Delta_cbf = 0, isolating the value of
%           the delay-robust tightening of Theorem 2.
%     'A4'  Full proposed DSC-SEC-HOCBF.
%     'A5'  Event-triggered PI-AGC comparator: transmission occurs when the
%           deviation from the last transmitted state exceeds a relative plus
%           absolute threshold, or when a maximum inter-event interval
%           elapses.  This is the standard event-triggered LFC baseline; the
%           maximum interval is what guarantees liveness.
%     'A6'  DELAY-COMPENSATED PREDICTIVE PI-AGC.  The central law is applied
%           not to the delayed state but to a forward prediction of it over
%           the measured age of data, obtained by integrating the nominal
%           model.  This is the conventional delay-compensation structure in
%           the load-frequency-control literature.
%
%           This arm exists because Reviewer 1 comment 11, Reviewer 2 comment
%           7 and Reviewer 3 comment 3 each ask for a comparison against a
%           delay-tolerant or predictive method, and an event-triggered
%           comparator alone does not answer them.  Reviewer 2 comment 10 asks
%           for the LIMITATION arising from such an absence to be stated; it
%           does not licence the absence, and treating it as though it did
%           would leave three comments unaddressed.
%
%           What remains outside scope is a full receding-horizon MPC with an
%           explicit cost and terminal set.  A6 is a one-step predictor, not
%           an optimisation over a horizon, and the paper must say so.  A full delay-tolerant MPC
%           comparator is NOT implemented and its absence is a stated
%           limitation, per Reviewer 2 comment 10; claiming otherwise would be
%           worse than admitting the gap.
%
%   INPUTS
%     arm    one of the labels above
%     t      current time [s]
%     x      TRUE local state, fresh at the edge every controller period
%     x_wan  state as delivered over the WAN, i.e. delayed; used by the
%            central AGC and by arm A2
%     rocof  measured d(df)/dt [p.u./s]
%     net    struct from the cyber layer: .alpha .pkt .age
%     p      parameters
%     ctx    .a_d .a_d_pred .h_min .h_max .cfg .Delta_cbf .qp (options)
%     cs     controller state, or [] to initialise
%
%   OUTPUTS
%     ug      central AGC command to the governor (zero when isolated)
%     u_edge  2x1 edge command actually applied
%     cs      updated state
%     ci      diagnostics, including .sigma and the QP record when applicable
%
%   The switching law applied here is the corrected one,
%       u_ref = alpha*(1-sigma)*u_AGC + (1 - alpha*(1-sigma))*u_loc
%   verified as E01 in verify_etm.  The submitted Eq. (36) used alpha*sigma
%   and selected the local law during healthy operation.  Note that the law
%   now selects the REFERENCE; the safety filter is applied unconditionally
%   on top of it.
%
%   THE ARMS ARE NOT A STRICTLY NESTED CHAIN, so a single linear attribution
%   of the A4-A0 difference is not meaningful.  They are a factor design:
%     arm  reference source        filter  tightening
%     A0   delayed AGC, no edge      no        -
%     A1   local PI                  no        -
%     A2   stale AGC                 yes       yes
%     A3   DA-ETM switched           yes       no
%     A4   DA-ETM switched           yes       yes
%     A5   event-triggered AGC       no        -
%   Report pairwise differences that isolate one factor, not a chain sum.
%
%   See also CTRL_HOCBF_QP, CTRL_PI_AGC, DA_ETM, VERIFY_QP.

if isempty(cs)
    cs = struct('qp',[],'etm',[],'u_prev',[0;0], ...
                'x_last_tx',x,'n_tx',0,'sigma',0, ...
                'lam_prev',[],'u_ref_prev',[0;0],'bump',[0;0],'t_bump',-inf, ...
                'ug_prev',0);
end
ci = struct('arm',arm,'sigma',0,'qi',[],'mech','none','T_pred',NaN,'tau_pred',NaN);

switch upper(arm)

%% ------------------------------------------------------------------ %%
case 'A0'                                   % legacy baseline, no edge
    [ug, ~] = ctrl_pi_agc(x_wan, p, 'agc_only');
    ug     = net.alpha*ug;
    u_edge = [0;0];

%% ------------------------------------------------------------------ %%
case 'A1'                                   % baseline WITH equal authority
    %  The central AGC runs on the DELAYED WAN state, as it must.  The edge
    %  dispatch runs on the LOCAL state, which is fresh at every controller
    %  period and survives a WAN outage.
    %
    %  Driving the edge law from x_wan instead makes the actuators sit idle
    %  during the attack, because the stale state shows ACE near zero.  A1
    %  would then be indistinguishable from A0 and the ablation would credit
    %  the switching with a gain that is really the actuators, which is the
    %  opposite of what this arm exists to measure.
    ug     = net.alpha*ctrl_pi_agc(x_wan, p, 'agc_only');
    [~,ue] = ctrl_pi_agc(x, p, 'edge_only');
    u_edge = ue;

%% ------------------------------------------------------------------ %%
case 'A5'                                   % event-triggered PI-AGC
    %  Relative deviation trigger with an absolute floor and a MAXIMUM
    %  inter-transmission interval.
    %
    %  A purely relative test dev > sigma*||x|| does not settle here.  The
    %  plant integrator z grows while the held state is stale, so ||x|| grows,
    %  the threshold grows with it, and the trigger stops firing: the arm then
    %  tracks a permanently outdated state and its frequency never returns to
    %  nominal.  A maximum inter-event time is standard in event-triggered
    %  control precisely to guarantee this liveness, and the deviation is
    %  measured against the TRANSMITTED state rather than the current one,
    %  which is the conventional form.
    sig_e = getfielddef(ctx,'etc_sigma',0.05);
    abs_e = getfielddef(ctx,'etc_abs',  1e-3);
    Tmax_e= getfielddef(ctx,'etc_Tmax', 1.0);
    if ~isfield(cs,'t_last_tx') || isempty(cs.t_last_tx), cs.t_last_tx = -inf; end

    %  The deviation is measured on the OUTPUT the scheme actually regulates,
    %  not on the full state.  The state vector contains the plant integrator
    %  z, which grows without bound while the held state is stale; including
    %  it makes the relative threshold enormous, the deviation test never
    %  fires, and only the maximum inter-event timer transmits.  The result is
    %  a sustained sawtooth at the timer period rather than a settled
    %  response, which is what a full-state trigger produced here.  Triggering
    %  on the measured output is also the conventional form in the
    %  event-triggered load-frequency-control literature.
    io   = [p.spec.ix.df, p.spec.ix.dPtie];
    dev  = norm(x(io) - cs.x_last_tx(io));
    thr  = sig_e*max(norm(cs.x_last_tx(io)), 1e-9) + abs_e;
    fire = (dev > thr) || ((t - cs.t_last_tx) >= Tmax_e);
    if fire && net.alpha
        cs.x_last_tx = x;
        cs.t_last_tx = t;
        cs.n_tx      = cs.n_tx + 1;
    end
    [ug, ue] = ctrl_pi_agc(cs.x_last_tx, p, 'agc_edge');
    ug     = net.alpha*ug;
    u_edge = ue;
    ci.n_tx = cs.n_tx;

%% ------------------------------------------------------------------ %%
case 'A6'                                   % delay-compensated predictive AGC
    %  The delayed state is propagated forward over the measured age of data
    %  before the central law is applied.  Prediction uses the NOMINAL model
    %  with the last applied commands held, which is what a delay compensator
    %  can actually know; it has no access to the true current state.
    %
    %  During an outage the age of data grows without bound, so the horizon is
    %  capped at the buffer limit.  Beyond that cap the predictor is
    %  extrapolating far outside any interval over which the frozen-input
    %  assumption holds, and reporting it as a fair comparator would flatter
    %  the proposed method rather than test it.
    tau_p = min(net.age, p.cyber.tau_pred_max);
    xp    = x_wan;
    if tau_p > 0
        nP = max(1, ceil(tau_p/p.spec.num.Ts_plant));
        hP = tau_p/nP;
        dP_hat = [0; 0; 0];                 % nominal: no disturbance forecast
        for q = 1:nP
            xp = cpmg_integrate(xp, cs.ug_prev, cs.u_prev, dP_hat, p, 0, hP, 'rk4');
        end
    end
    ug     = net.alpha*ctrl_pi_agc(xp, p, 'agc_only');   % central path compensated
    [~,ue] = ctrl_pi_agc(x,  p, 'edge_only');            % edge on LOCAL state, as A1
    u_edge = ue;
    cs.ug_prev = ug;
    ci.tau_pred = tau_p;

%% ------------------------------------------------------------------ %%
case {'A2','A3','A4'}
    %  ARCHITECTURE.  The safety filter is applied at EVERY controller period,
    %  on the LOCAL state, which is fresh regardless of the WAN.  What the
    %  event trigger selects is the REFERENCE the filter tracks, not whether
    %  the filter runs.
    %
    %  An earlier version applied u_ref unfiltered whenever the link was
    %  healthy and used the QP output only while isolated.  That leaves the
    %  system with no safety guarantee during normal operation, and it made
    %  arm A2 change two things at once (barrier added AND the barrier fed a
    %  stale state), which is why the linear attribution came out as
    %  "barrier -0.68, switching +0.93": the same effect counted twice with
    %  opposite signs.  Feeding the barrier a delayed state was never
    %  physically motivated either, since an edge controller always has local
    %  measurements.
    %
    %  A2 has NO staleness detection, so it keeps tracking the last AGC
    %  command it received, however old.  A3 and A4 switch the reference to a
    %  locally computed one once the trigger fires.

    if strcmpi(arm,'A2')
        sigma = 0;  ci.mech = 'disabled';
    else
        ectx = struct('a_d',ctx.a_d,'a_d_pred',ctx.a_d_pred, ...
                      'h_min',ctx.h_min,'h_max',ctx.h_max,'cfg',ctx.cfg, ...
                      'tie_in',getfielddef(ctx,'tie_in',0));
        [sigma, cs.etm, edi] = da_etm(t, x, rocof, net.pkt, net.age, p, ectx, cs.etm);
        ci.mech   = edi.mech;
        ci.T_pred = edi.T_pred;
    end
    cs.sigma = sigma;  ci.sigma = sigma;

    % ---- reference selection, the corrected switching law ---------------
    %      u_ref = alpha*(1-sigma)*u_AGC + (1 - alpha*(1-sigma))*u_local
    %  The isolation flag is passed to both reference laws so that the charge
    %  restoration term is disabled whenever the edge is isolated, per the
    %  timescale-separation rule in ctrl_pi_agc.
    gs = struct('sigma', sigma);
    [ug_raw, u_agc] = ctrl_pi_agc(x_wan, p, 'agc_edge', gs);   % delayed, WAN side
    [~,      u_loc] = ctrl_pi_agc(x,     p, 'edge_only', gs);  % local, fresh

    if strcmpi(arm,'A2')
        lam = 1;                       % no detection: keep using the stale AGC
    else
        lam = net.alpha*(1 - sigma);
    end
    %  ---- who drives the GOVERNOR --------------------------------------
    %  The AGC reaches the governor whenever the LINK is up.  sigma decides
    %  who drives the EDGE ACTUATORS; it must not sever secondary control to
    %  the turbine on a healthy link.  Using alpha*(1-sigma) here meant a
    %  predictive trigger disconnected a perfectly good AGC, leaving the
    %  governor with droop only and forcing the storage to carry the entire
    %  steady-state imbalance until it drained.
    ug = net.alpha*ug_raw;

    %  ---- bumpless transfer on the edge reference -----------------------
    %  zACE is a PLANT state shared by both reference laws and neither can
    %  reset it.  During isolation the local law drives df to zero while zACE
    %  keeps whatever it accumulated, so switching back steps the reference,
    %  the step produces a RoCoF spike, and the spike re-fires the predictive
    %  trigger.  The result is a limit cycle: 58 isolation episodes over 180 s
    %  with the frequency oscillating indefinitely.
    %
    %  The offset below makes the reference continuous across every switch and
    %  then decays it with time constant T_bump, which is standard bumpless
    %  transfer between two controllers that cannot share an integrator.
    u_raw = lam*u_agc + (1 - lam)*u_loc;
    if isempty(cs.lam_prev) || abs(lam - cs.lam_prev) > 1e-9
        cs.bump   = cs.u_ref_prev - u_raw;
        cs.t_bump = t;
    end
    dec   = exp(-(t - cs.t_bump)/p.etm.T_bump);
    u_ref = u_raw + cs.bump*dec;

    cs.lam_prev   = lam;
    cs.u_ref_prev = u_ref;

    % ---- delay tightening ------------------------------------------------
    switch upper(arm)
        case 'A3', Delta = 0;               % ablate Theorem 2
        otherwise, Delta = ctx.Delta_cbf;
    end

    % ---- safety filter, ALWAYS applied, on the LOCAL state ---------------
    qopt = ctx.qp;  qopt.Delta_cbf = Delta;
    qopt.tie_in   = getfielddef(ctx,'tie_in',0);
    [u_edge, cs.qp, qi] = ctrl_hocbf_qp(x, rocof, u_ref, p, cs.qp, qopt);
    ci.qi    = qi;
    ci.u_ref = u_ref;
    ci.lam   = lam;

otherwise
    error('ctrl_bank:arm','Unknown ablation arm "%s"', arm);
end

u_edge = [ min(max(u_edge(1), p.uB_min), p.uB_max) ;
           min(max(u_edge(2), p.uA_min), p.uA_max) ];
cs.u_prev = u_edge;
ci.u_edge = u_edge;
ci.ug     = ug;
end

%% ==================================================================== %%
function v = getfielddef(s, f, dflt)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

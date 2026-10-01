function [sigma, es, di] = da_etm(t, x, rocof, pkt_arrived, age, p, ctx, es)
%DA_ETM  Delay-Aware Event-Triggering Mechanism.
%
%   [sigma, es, di] = DA_ETM(t, x, rocof, pkt_arrived, age, p, ctx, es)
%
%   sigma = 0  trust the WAN, apply the central AGC command
%   sigma = 1  isolated, apply the local barrier-filtered command
%
%   WHY A BARE AGE-OF-DATA THRESHOLD IS NOT ENOUGH
%   ----------------------------------------------
%   Layer 4 measured the time from DoS onset to loss of the barrier under a
%   frozen edge command.  For small initial margins that time is shorter than
%   any achievable packet-loss detection latency: a link with period T_pkt
%   cannot declare a loss faster than n_miss*T_pkt, yet the barrier can be lost
%   sooner than that.  A trigger that waits for evidence of an attack is
%   therefore structurally too slow for the states that need it most.
%
%   The mechanism implemented here has three paths, and the ordering matters.
%
%     M1  PREDICTIVE (safety-critical).  Using only LOCAL measurements, which
%         are fresh at every controller period regardless of the WAN, predict
%         the time for the barrier to reach the feasibility floor UNDER THE
%         HYPOTHESIS THAT THE EDGE COMMAND FREEZES NOW.  If that is shorter
%         than the reaction horizon, isolate.  This path does not wait for the
%         attack to be detected; it responds to the margin.
%
%         Engaging early is not a false positive.  Under sigma = 1 the edge
%         runs the HO-CBF-QP with the AGC command as its reference, so no
%         function is lost; only the unfiltered pass-through is given up.
%
%     M2  PACKET LOSS.  n_miss consecutive missed arrivals declare a DoS.
%         This is the path that actually identifies the attack, and it is
%         deterministic, but it is bounded below by the link period.
%
%     M3  AGE-OF-DATA BACKSTOP.  age > tau_trigger catches slow degradation
%         that neither of the above sees.  tau_trigger sits strictly between
%         the worst healthy age and the delay feasibility limit tau_max.
%
%   RECONNECTION requires all of: the dwell time elapsed, n_ok consecutive
%   fresh packets, and the margin recovered above a stated fraction of the
%   safe set.  Without the margin condition the mechanism chatters at the
%   barrier, switching authority away from the local controller exactly when
%   it is needed.
%
%   INPUTS
%     t            current time [s]
%     x            LOCAL state, fresh at every controller period
%     rocof        measured d(df)/dt [p.u./s].  REQUIRED.  The barrier
%                  derivative depends on the unmeasured disturbance, so
%                  evaluating it with a zero disturbance reports an
%                  apparently healthy margin rate while the frequency is
%                  collapsing. See estimate_disturbance.
%     pkt_arrived  true when a WAN packet arrived at this instant
%     age          current age of data [s]
%     p            parameters
%     ctx          struct: .a_d      curvature bound from Layer 4
%                          .h_min    feasibility floor from Layer 4
%                          .h_max    -df_min, the full safe set width
%                          .cfg      barrier class-K config
%     es           mechanism state, or [] to initialise
%
%   OUTPUTS
%     sigma  0 or 1
%     es     updated state
%     di     diagnostics: .mech which path fired, .T_pred, .h, .hdot,
%            .n_missed, .dwell_ok, .n_fresh
%
%   See also TIME_TO_FLOOR, DOS_GENERATOR, DELAY_BUFFER, VERIFY_ETM.

if isempty(es)
    es = struct('sigma',0,'t_iso',-inf,'t_last_pkt',t,'n_fresh',0, ...
                'n_trig',0,'mech',{{}},'t_trig',[],'mech_last','none');
end

ctx = setdef(ctx,'a_d',   1.7960);
%   TWO DIFFERENT CURVATURES, AND USING THE WRONG ONE BREAKS THE MECHANISM.
%
%   ctx.a_d       ADVERSARIAL bound: sup over the operating domain of -d2h/dt2
%                 with the control input ranging over its whole box.  This is
%                 the right quantity for a SAFETY BOUND, and compute_scatd
%                 uses it.  It is the WRONG quantity here.
%
%   ctx.a_d_pred  FROZEN-INPUT curvature, measured on the trajectories the
%                 trigger is actually predicting.  The trigger asks "if the
%                 WAN is lost now and the edge command freezes, how long do I
%                 have", and that hypothetical has no adversary in it.
%
%   With the adversarial value the predictor is pessimistic everywhere: it
%   fires under benign healthy transients, and the reconnection condition
%   T_pred >= 3*T_react becomes unsatisfiable at ANY state, so the mechanism
%   latches isolated forever.  Both were observed.
%
%   a_d_pred is a TRIGGER DESIGN parameter validated by verify_etm E05.  It is
%   NOT a safety bound and must never be used to support a safety claim; the
%   invariance guarantee comes from the QP, not from this prediction.
ctx = setdef(ctx,'a_d_pred', ctx.a_d);
ctx = setdef(ctx,'h_min', 0);
ctx = setdef(ctx,'h_max', -p.spec.safety.df_min_pu);
ctx = setdef(ctx,'cfg',   struct('a_f',[8 8],'Delta_cbf',0));

%% ------------------------------------------------------------------ %%
%  Local barrier state, always fresh
%% ------------------------------------------------------------------ %%
%   The disturbance is reconstructed algebraically from measured RoCoF, so
%   that Lf h below is the TRUE margin rate rather than a disturbance-free
%   estimate of it.  Passing zeros here is the defect that gave the predictive
%   trigger a timing margin of exactly zero.
d_hat = estimate_disturbance(x, rocof, p);
B     = barrier_functions(x, 0, d_hat, p, 0, ctx.cfg);
h     = B.f.h;
hdot  = B.f.Lfh;

T_pred = time_to_floor(h, hdot, ctx.a_d_pred, ctx.h_min);

%% ------------------------------------------------------------------ %%
%  Packet bookkeeping
%% ------------------------------------------------------------------ %%
%   n_fresh counts consecutive received PACKETS, not consecutive controller
%   steps on which a packet happened to land.  da_etm is called every
%   controller period (Ts_ctrl) while packets arrive every T_pkt, so with
%   T_pkt = 50 ms and Ts_ctrl = 10 ms an arrival occurs on one call in five.
%   Resetting the counter on every call without an arrival caps it at 1,
%   which makes fresh_ok unsatisfiable for any n_ok > 1 and latches the
%   mechanism isolated permanently.  The counter is reset only when a packet
%   is actually MISSED.
if pkt_arrived
    es.t_last_pkt = t;
    es.n_fresh    = es.n_fresh + 1;
elseif (t - es.t_last_pkt) > 1.5*p.cyber.T_pkt
    es.n_fresh    = 0;                 % a packet was missed, not merely absent
end
n_missed = floor((t - es.t_last_pkt)/p.cyber.T_pkt);

%% ------------------------------------------------------------------ %%
%  Trigger evaluation
%% ------------------------------------------------------------------ %%
m1 = T_pred   <  p.etm.T_react;          % predictive, safety-critical
m2 = n_missed >= p.etm.n_miss;           % packet loss
m3 = age      >  p.etm.tau_trigger;      % age-of-data backstop

mech = 'none';
if es.sigma == 0
    if m1
        mech = 'predictive';
    elseif m2
        mech = 'packet_loss';
    elseif m3
        mech = 'age_of_data';
    end
    if ~strcmp(mech,'none')
        es.sigma     = 1;
        es.t_iso     = t;
        es.n_trig    = es.n_trig + 1;
        es.mech{end+1} = mech;
        es.t_trig(end+1) = t;
        es.mech_last = mech;
    end
else
    %% -------------------------------------------------------------- %%
    %  Reconnection.  All conditions, no exceptions.
    %% -------------------------------------------------------------- %%
    dwell_ok  = (t - es.t_iso) >= p.etm.T_dwell;
    fresh_ok  = es.n_fresh     >= p.etm.n_ok;
    hi_ok     = h      >= p.etm.h_recover*ctx.h_max;
    Tp_ok     = T_pred >= p.etm.k_recover*p.etm.T_react;
    age_ok    = age    <= p.etm.tau_trigger;
    margin_ok = hi_ok && Tp_ok && age_ok;
    if dwell_ok && fresh_ok && margin_ok
        es.sigma = 0;
        mech = 'reconnect';
        es.mech_last = mech;
    end
end

sigma = es.sigma;

%   EVERY reconnection sub-condition is returned, so that a mechanism which
%   fails to release can be diagnosed from a log instead of by guesswork.
%   Three separate latch defects have now been traced to a single unreachable
%   condition, and in each case the blocking term was identified only after a
%   trace was instrumented.
di = struct('mech',mech,'T_pred',T_pred,'h',h,'hdot',hdot,'w_hat',d_hat(1), ...
            'n_missed',n_missed,'n_fresh',es.n_fresh, ...
            'age',age,'m1',m1,'m2',m2,'m3',m3, ...
            'dwell',t - es.t_iso, ...
            'rc_dwell', (t - es.t_iso) >= p.etm.T_dwell, ...
            'rc_fresh', es.n_fresh >= p.etm.n_ok, ...
            'rc_h',     h      >= p.etm.h_recover*ctx.h_max, ...
            'rc_Tpred', T_pred >= p.etm.k_recover*p.etm.T_react, ...
            'rc_age',   age    <= p.etm.tau_trigger);
end

%% ==================================================================== %%
function s = setdef(s,f,v)
if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end

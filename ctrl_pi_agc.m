function [ug, u_edge, ci] = ctrl_pi_agc(x, p, mode, gains)
%CTRL_PI_AGC  Conventional PI area control error controller.
%
%   [ug, u_edge, ci] = CTRL_PI_AGC(x, p, mode, gains)
%
%   Serves three purposes in this repository:
%     1. it is the nominal reference law u_ref consumed by the HO-CBF-QP;
%     2. it is ablation arm A0, the legacy baseline with NO edge actuation;
%     3. it is ablation arm A1, the SAME baseline given IDENTICAL BESS and
%        AIAC authority to the proposed controller.
%
%   ARM A1 IS THE ONE THAT MATTERS.  The previous submission compared the
%   proposed framework only against a PI-AGC that had no BESS and no AIAC at
%   all, so none of the reported improvement could be attributed to the
%   controller rather than to the extra actuators.  That is Reviewer 3
%   comment 3 and Reviewer 1 comment 11.  Arm A1 isolates the question.
%
%   INPUTS
%     x      9x1     canonical state
%     p      struct  parameters
%     mode   char    'agc_only'  -> A0: central AGC only, u_edge = 0
%                    'agc_edge'  -> A1: central AGC plus PI edge dispatch
%                    'edge_only' -> local PI edge law with no central AGC,
%                                   used as the u_ref supplied to the QP
%     gains  struct  optional overrides: .Kp .Ki .Kp_e .Ki_e
%
%   OUTPUTS
%     ug      1x1    central AGC command to the governor reference    [p.u.]
%     u_edge  2x1    [u_B ; u_A] edge commands                        [p.u.]
%     ci      struct .ACE .zACE .sat_g .sat_B .sat_A saturation flags
%
%   SIGN CONVENTION (D3)
%     ACE = beta*df + dPtie.  Under-frequency gives ACE < 0.
%     u_B = -(Kp*ACE + Ki*zACE) > 0  : BESS discharges, injecting power.
%     u_A = +(Kp*ACE + Ki*zACE) < 0  : AIAC sheds cooling load.
%     The opposite signs are correct and follow from the AIAC entering the
%     swing equation as a consumption term.
%
%   INTEGRAL WIND-UP
%     zACE is a PLANT state, not a controller state, so the controller cannot
%     reset it.  Output saturation is applied and the saturation flags are
%     returned so that wind-up episodes can be reported rather than hidden.
%     The edge channels carry no integral action at all (see below), which
%     also removes them from this concern.
%
%   See also BARRIER_FUNCTIONS, ESTIMATE_CONTRACTION_RATES, CTRL_BANK.

if nargin < 3 || isempty(mode),  mode  = 'agc_only'; end
if nargin < 4,                   gains = struct();   end

g = gains;
if ~isfield(g,'Kp')   || isempty(g.Kp),   g.Kp   = p.agc.Kp;  end
if ~isfield(g,'Ki')   || isempty(g.Ki),   g.Ki   = p.agc.Ki;  end
if ~isfield(g,'Kp_e') || isempty(g.Kp_e), g.Kp_e = p.agc.Kp_e;  end
if ~isfield(g,'Ki_e') || isempty(g.Ki_e), g.Ki_e = p.agc.Ki_e;  end
if ~isfield(g,'K_soc')|| isempty(g.K_soc),g.K_soc= p.agc.K_soc; end
%  sigma is the isolation flag of the event-triggering mechanism.  Arms with
%  no trigger pass zero, which is correct: they are never isolated.
if ~isfield(g,'sigma')|| isempty(g.sigma),g.sigma = 0;          end

ix   = p.spec.ix;
df   = x(ix.df);
dPti = x(ix.dPtie);
zACE = x(ix.zACE);

ACE = p.beta*df + dPti;

%% ------------------------------------------------------------------ %%
%  Central AGC to the governor reference
%% ------------------------------------------------------------------ %%
ug_raw = -( g.Kp*ACE + g.Ki*zACE );
ug     = min(max(ug_raw, p.agc.u_min), p.agc.u_max);
sat_g  = abs(ug - ug_raw) > eps;

%% ------------------------------------------------------------------ %%
%  Edge dispatch
%% ------------------------------------------------------------------ %%
%  NO INTEGRAL ACTION ON THE EDGE REFERENCE, and a state-of-charge
%  restoration term on the storage channel.
%
%  zACE is the PLANT integral of ACE.  After a load step ACE returns to zero
%  but zACE settles at a nonzero constant, because that constant is what
%  holds the governor's steady-state output.  An edge law of the form
%  -(Kp_e*ACE + Ki_e*zACE) therefore keeps commanding -Ki_e*zACE forever, so
%  the storage is asked to supply steady-state ENERGY: in a 180 s run after a
%  2 s attack the battery discharged continuously from SoC 0.50 to 0.36 with
%  the link healthy and the frequency exactly at nominal.
%
%  The correct division of labour is that the governor and the AGC hold the
%  steady state while the storage handles the transient and then recovers.
%  Ki_e therefore defaults to zero and a restoration term returns the BESS to
%  its nominal charge.  The barrier still commands whatever the transient
%  needs; this law only sets the nominal reference the filter tracks.
%
%  No restoration term is applied to the AIAC channel: its state of charge
%  analogue is the building temperature, which has its own barrier.
SoC = x(ix.dSB);
%  TIMESCALE SEPARATION OF THE CHARGE RESTORATION TERM.
%
%  Unmodulated, the restoration term biases the reference in BOTH directions
%  depending on the initial charge, and neither direction is admissible during
%  a frequency excursion.  Above nominal charge it commands charging while the
%  frequency falls, so a higher starting charge yields a SMALLER barrier
%  margin.  Below nominal it commands discharge, which flatters the margin for
%  the same wrong reason.  Every result taken at a charge away from nominal is
%  therefore biased in one direction or the other, the low-charge points
%  optimistically and the high-charge points pessimistically.
%
%  State-of-charge restoration is a tertiary-timescale energy-management
%  function and must not act on the primary or secondary reference while the
%  frequency is disturbed.  The term is therefore enabled only when the
%  frequency is inside the governor dead-band AND the edge is not isolated:
%
%      u_B^ref = -Kp_e*ACE - K_soc*(S - S_nom) * 1{ |df| <= df_db  and  sigma = 0 }
%
%  Both gating signals already exist in the model, so this introduces no new
%  parameter and no new tuning. It is the same class of correction as setting
%  Ki_e to zero, and it rests on timescale separation rather than on a fitted
%  threshold. Phi(SoC) continues to enforce the charge limits, so nothing in
%  the barrier construction changes.
quiet  = abs(x(ix.df)) <= p.df_db;
g_soc  = double(quiet && (g.sigma == 0));

uB_raw = -( g.Kp_e*ACE + g.Ki_e*zACE ) - g_soc*g.K_soc*(SoC - p.SoC_nom);
uA_raw = +( g.Kp_e*ACE + g.Ki_e*zACE );

uB = min(max(uB_raw, p.uB_min), p.uB_max);
uA = min(max(uA_raw, p.uA_min), p.uA_max);

switch lower(mode)
    case 'agc_only'                     % A0
        u_edge = [0; 0];
        sat_B  = false;  sat_A = false;
    case 'agc_edge'                     % A1
        u_edge = [uB; uA];
        sat_B  = abs(uB - uB_raw) > eps;
        sat_A  = abs(uA - uA_raw) > eps;
    case 'edge_only'                    % u_ref for the QP
        ug     = 0;  sat_g = false;
        u_edge = [uB; uA];
        sat_B  = abs(uB - uB_raw) > eps;
        sat_A  = abs(uA - uA_raw) > eps;
    otherwise
        error('ctrl_pi_agc:mode','Unknown mode "%s"', mode);
end

ci = struct('ACE',ACE,'zACE',zACE,'SoC_err',SoC - p.SoC_nom,'sat_g',sat_g,'sat_B',sat_B,'sat_A',sat_A, ...
            'ug_raw',ug_raw,'uB_raw',uB_raw,'uA_raw',uA_raw,'mode',mode, ...
            'g_soc',g_soc);
end

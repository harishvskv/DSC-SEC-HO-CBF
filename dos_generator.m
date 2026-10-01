function A = dos_generator(p, t, mode, opt)
%DOS_GENERATOR  DoS attack sequence and WAN latency profile.
%
%   A = DOS_GENERATOR(p, t, mode, opt) returns the network availability signal,
%   the packet arrival schedule and the latency profile over the time grid t.
%
%   mode = 'none'          healthy network throughout
%          'deterministic' attacks at specified start times and durations
%          'stochastic'    seeded random attacks, for the Monte Carlo study
%
%   EVERY DRAW IS RECORDED.  A.draws holds the seed, the realised attack count,
%   start times and durations, and the distribution family and support of each
%   random quantity.  Reviewer 1 comment 6 asks for exactly this, and it cannot
%   be supplied after the fact.
%
%   DISTRIBUTIONS ARE TRUNCATED AND SAID TO BE
%     The legacy study drew inertia, damping and time constants from
%     UNBOUNDED Gaussians while the manuscript described them as plus or minus
%     a percentage, which reads as a uniform support.  Roughly five per cent of
%     those draws lay outside the stated range and a small tail admitted
%     non-physical values.  Every distribution here is either uniform on a
%     stated interval or a Gaussian truncated to a stated interval, and the
%     choice is recorded in A.draws.spec.
%
%   OUTPUTS (struct A)
%     A.alpha       1xN logical, true when the WAN is available
%     A.tau_net     1xN WAN latency [s]
%     A.pkt         1xN logical, true at instants when a packet ARRIVES
%     A.t_stamp     1xN timestamp carried by the most recent arrival
%     A.age         1xN age of data seen by the edge = t - t_stamp
%     A.windows     nA x 2 attack [start stop]
%     A.draws       provenance record
%
%   TIMING CONSISTENCY
%     The worst healthy age of data is p.cyber.T_pkt + p.cyber.tau_nom_hi, and
%     the constructor asserts that it lies strictly below p.etm.tau_trigger,
%     which in turn must lie below the delay feasibility limit tau_max.  If
%     that ordering is violated the trigger either never fires or fires
%     constantly on healthy jitter.
%
%   See also DA_ETM, DELAY_BUFFER, VERIFY_ETM.

if nargin < 3 || isempty(mode), mode = 'none'; end
if nargin < 4, opt = struct(); end
opt = setdef(opt,'seed',1);
opt = setdef(opt,'starts',120);
opt = setdef(opt,'durations',3.0);
opt = setdef(opt,'n_range',[1 3]);
opt = setdef(opt,'start_range',[115 125]);
opt = setdef(opt,'dur_range',[2.0 3.4]);
opt = setdef(opt,'spike_range',[0.5 2.0]);
opt = setdef(opt,'jitter','uniform');

t = t(:)';
N = numel(t);
dt = t(2) - t(1);

%% ------------------------------------------------------------------ %%
%  Timing consistency assertions
%% ------------------------------------------------------------------ %%
age_healthy = p.cyber.T_pkt + p.cyber.tau_nom_hi;
assert(age_healthy < p.etm.tau_trigger, ...
    ['dos_generator: worst healthy age of data %.4f s is not below the ' ...
     'event trigger %.4f s. The trigger would fire on normal jitter.'], ...
     age_healthy, p.etm.tau_trigger);

%% ------------------------------------------------------------------ %%
%  Attack windows
%% ------------------------------------------------------------------ %%
draws = struct('mode',mode,'seed',opt.seed);
switch lower(mode)
    case 'none'
        W = zeros(0,2);
        draws.spec = 'no attack';

    case 'deterministic'
        s0 = opt.starts(:);  d0 = opt.durations(:);
        if isscalar(d0), d0 = repmat(d0,size(s0)); end
        W = [s0, s0 + d0];
        draws.spec = 'deterministic start times and durations';

    case 'stochastic'
        rng(opt.seed,'twister');
        nA = randi(opt.n_range);
        s0 = opt.start_range(1) + diff(opt.start_range)*rand(nA,1);
        d0 = opt.dur_range(1)   + diff(opt.dur_range)*rand(nA,1);
        W  = [s0, s0 + d0];
        draws.spec = struct( ...
            'count',    sprintf('uniform integer on [%d, %d]', opt.n_range), ...
            'start',    sprintf('uniform on [%.2f, %.2f] s', opt.start_range), ...
            'duration', sprintf('uniform on [%.2f, %.2f] s', opt.dur_range), ...
            'spike',    sprintf('uniform on [%.2f, %.2f] s', opt.spike_range));
        draws.count = nA;

    otherwise
        error('dos_generator:mode','Unknown mode "%s"', mode);
end
draws.windows = W;

alpha = true(1,N);
for k = 1:size(W,1)
    alpha(t >= W(k,1) & t <= W(k,2)) = false;
end

%% ------------------------------------------------------------------ %%
%  WAN latency
%% ------------------------------------------------------------------ %%
switch lower(mode)
    case 'stochastic'
        tau_net = p.cyber.tau_nom_lo + ...
                  (p.cyber.tau_nom_hi - p.cyber.tau_nom_lo)*rand(1,N);
        spike = opt.spike_range(1) + diff(opt.spike_range)*rand(1,N);
        tau_net(~alpha) = tau_net(~alpha) + spike(~alpha);
    otherwise
        % Deterministic sinusoidal jitter inside the declared band
        mid = 0.5*(p.cyber.tau_nom_lo + p.cyber.tau_nom_hi);
        amp = 0.5*(p.cyber.tau_nom_hi - p.cyber.tau_nom_lo);
        tau_net = mid + amp*sin(0.5*t);
end
% Cap the latency at the declared spike support, NOT at the age of data.
% Those are different quantities and conflating them silently truncated the
% spike in the legacy study.
tau_net = min(tau_net, opt.spike_range(2) + p.cyber.tau_nom_hi);

%% ------------------------------------------------------------------ %%
%  Packet arrivals and age of data
%% ------------------------------------------------------------------ %%
%   A packet is generated every T_pkt.  It ARRIVES tau_net later, and only if
%   the WAN was available over the whole flight.  The timestamp it carries is
%   the generation instant, so the age of data at the edge is t - t_gen.
pkt     = false(1,N);
t_stamp = nan(1,N);
kGen    = round(p.cyber.T_pkt/dt);
kGen    = max(kGen,1);

last_stamp = t(1);
for n = 1:N
    if mod(n-1, kGen) == 0
        tg = t(n);
        na = n + round(tau_net(n)/dt);
        if na <= N && all(alpha(n:min(na,N)))
            pkt(na) = true;
            t_stamp(na) = tg;
        end
    end
end
for n = 1:N
    if pkt(n), last_stamp = t_stamp(n); end
    t_stamp(n) = last_stamp;
end
age = t - t_stamp;

%% ------------------------------------------------------------------ %%
A = struct('t',t,'alpha',alpha,'tau_net',tau_net,'pkt',pkt, ...
           't_stamp',t_stamp,'age',age,'windows',W,'draws',draws, ...
           'age_healthy_max',age_healthy,'T_pkt',p.cyber.T_pkt);

A.age_max_observed = max(age);
% The age of data produced by an attack of duration D is
%     D + T_pkt + tau_nom_hi
% not D.  If this fires, raise p.cyber.dos_dur_max, which resizes both
% age_dos_max and tau_buffer consistently; do not patch tau_buffer alone.
A.dur_max_observed = 0;
if ~isempty(W), A.dur_max_observed = max(W(:,2) - W(:,1)); end
assert(A.dur_max_observed <= p.cyber.dos_dur_max + 1e-9, ...
    ['dos_generator: attack duration %.3f s exceeds p.cyber.dos_dur_max ' ...
     '= %.3f s. Raise dos_dur_max; the buffer horizon follows from it.'], ...
     A.dur_max_observed, p.cyber.dos_dur_max);
assert(A.age_max_observed <= p.cyber.tau_buffer, ...
    ['dos_generator: observed age of data %.3f s exceeds the delay buffer ' ...
     'horizon %.3f s (attack duration %.3f s + T_pkt %.3f + latency %.3f). ' ...
     'Raise p.cyber.dos_dur_max.'], ...
     A.age_max_observed, p.cyber.tau_buffer, A.dur_max_observed, ...
     p.cyber.T_pkt, p.cyber.tau_nom_hi);
end

function s = setdef(s,f,v)
if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end

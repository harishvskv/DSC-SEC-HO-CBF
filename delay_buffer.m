function varargout = delay_buffer(cmd, varargin)
%DELAY_BUFFER  Ring buffer of past states, for sampled-data and delayed access.
%
%   buf          = DELAY_BUFFER('init',  p, x0)
%   buf          = DELAY_BUFFER('push',  buf, t, x)
%   [x,age,info] = DELAY_BUFFER('fetch', buf, t_now, age_request)
%   [x,age,info] = DELAY_BUFFER('fetch_at', buf, t_stamp)
%   s            = DELAY_BUFFER('stats', buf)
%
%   WHY THIS FILE FAILS LOUDLY
%   --------------------------
%   The legacy implementation held 20 samples, 0.19 s of history, while the
%   stochastic study requested delays of up to 2.0 s.  The index arithmetic
%   was clamped with max(1, ...), so every delay above 0.19 s was SILENTLY
%   truncated to 0.19 s.  The DoS-induced delay spike, one of the headline
%   features of that study, never entered the dynamics at all and nothing
%   reported it.
%
%   Here the horizon is sized from p.cyber.tau_buffer, which params asserts is
%   larger than p.cyber.tau_dos_hi, and any request beyond the stored history
%   raises an error rather than returning the oldest sample.  A saturating
%   variant is available only by passing info.allow_clamp explicitly, and it
%   still sets info.clamped so the caller cannot ignore it.
%
%   FIELDS OF buf
%     .X        nx x N   circular store of states
%     .T        1  x N   matching timestamps
%     .head     index of the most recent entry
%     .count    number of valid entries
%     .N        capacity
%     .Ts       nominal push interval, used for index arithmetic
%     .horizon  Ts*(N-1), the oldest retrievable age
%     .maxAgeSeen, .nClamped  diagnostics
%
%   See also DA_ETM, DOS_GENERATOR, CPMG_INTEGRATE.

switch lower(cmd)

%% ------------------------------------------------------------------ %%
case 'init'
    p  = varargin{1};
    x0 = varargin{2};
    Ts = p.spec.num.Ts_plant;
    N  = ceil(p.cyber.tau_buffer/Ts) + 2;

    assert(p.cyber.tau_buffer > p.cyber.tau_dos_hi, ...
        ['delay_buffer: horizon %.3f s does not exceed the worst-case age of ' ...
         'data %.3f s. This is the exact condition the legacy code violated.'], ...
         p.cyber.tau_buffer, p.cyber.tau_dos_hi);

    buf.X = repmat(x0(:), 1, N);
    buf.T = -inf(1, N);
    buf.T(1) = 0;
    buf.head = 1;
    buf.count = 1;
    buf.N = N;
    buf.Ts = Ts;
    buf.horizon = Ts*(N-1);
    buf.maxAgeSeen = 0;
    buf.nClamped = 0;
    buf.nx = numel(x0);
    varargout{1} = buf;

%% ------------------------------------------------------------------ %%
case 'push'
    buf = varargin{1};
    t   = varargin{2};
    x   = varargin{3};
    buf.head = mod(buf.head, buf.N) + 1;
    buf.X(:, buf.head) = x(:);
    buf.T(buf.head) = t;
    buf.count = min(buf.count + 1, buf.N);
    varargout{1} = buf;

%% ------------------------------------------------------------------ %%
case 'fetch'
    buf   = varargin{1};
    t_now = varargin{2};
    age   = varargin{3};
    if nargin >= 5, allow_clamp = varargin{4}; else, allow_clamp = false; end

    buf.maxAgeSeen = max(buf.maxAgeSeen, age);

    k = round(age/buf.Ts);
    clamped = false;
    if k > buf.count - 1
        if allow_clamp
            k = buf.count - 1;
            clamped = true;
            buf.nClamped = buf.nClamped + 1;
        else
            error('delay_buffer:horizon', ...
                ['Requested age %.4f s exceeds the stored history %.4f s ' ...
                 '(capacity %d samples at Ts = %.4f s). Increase ' ...
                 'p.cyber.tau_buffer. Silently clamping here is the defect ' ...
                 'that voided the legacy stochastic study.'], ...
                 age, (buf.count-1)*buf.Ts, buf.N, buf.Ts);
        end
    end
    k = max(k, 0);

    idx = mod(buf.head - k - 1, buf.N) + 1;
    varargout{1} = buf.X(:, idx);
    varargout{2} = t_now - buf.T(idx);
    varargout{3} = struct('index',idx,'k',k,'clamped',clamped, ...
                          'stamp',buf.T(idx),'buf',buf);

%% ------------------------------------------------------------------ %%
case 'fetch_at'
    buf    = varargin{1};
    tstamp = varargin{2};
    valid  = isfinite(buf.T);
    [~, j] = min(abs(buf.T - tstamp) + 1e9*(~valid));
    varargout{1} = buf.X(:, j);
    varargout{2} = buf.T(buf.head) - buf.T(j);
    varargout{3} = struct('index',j,'k',NaN,'clamped',false, ...
                          'stamp',buf.T(j),'buf',buf);

%% ------------------------------------------------------------------ %%
case 'stats'
    buf = varargin{1};
    varargout{1} = struct('capacity',buf.N,'count',buf.count, ...
        'horizon',buf.horizon,'maxAgeSeen',buf.maxAgeSeen, ...
        'nClamped',buf.nClamped,'Ts',buf.Ts);

otherwise
    error('delay_buffer:cmd','Unknown command "%s"', cmd);
end
end

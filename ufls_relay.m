function [shed, st] = ufls_relay(df, t, st, spec)
%UFLS_RELAY  Staged under-frequency load shedding relay, IEGC 2023 settings.
%
%   [shed, st] = UFLS_RELAY(df, t, st, spec) models the four default UFR
%   stages (49.40 / 49.20 / 49.00 / 48.80 Hz) as a discrete element that is
%   EXTERNAL to the plant.  Keeping the relay outside cpmg_dynamics matters:
%   the legacy implementation buried load shedding inside the drift field as
%   a boolean flag, which made the vector field discontinuous in a way the
%   barrier module could not see.
%
%   INPUTS
%     df    1x1  frequency deviation                                 [p.u.]
%     t     1x1  current simulation time                             [s]
%     st    1x1  relay state struct, or [] to initialise
%     spec  1x1  cpmg_spec()
%
%   OUTPUTS
%     shed  1x1  cumulative fraction of load shed, in [0,1]
%     st    1x1  updated relay state
%
%   RELAY LOGIC
%     A stage arms when df falls below its pick-up level.  It trips after the
%     stage operating time has elapsed with df still below pick-up.  Once
%     tripped a stage latches; it resets only if frequency recovers above
%     spec.ufls.reset_Hz, and reset is optional (default: latched).
%
%   REPORTING
%     st.trip_time(k) is NaN until stage k trips.  st.armed_time records the
%     first arming instant of each stage, which is what the manuscript should
%     report as a near-miss indicator alongside the nadir.
%
%   See also CPMG_SPEC.

if nargin < 4 || isempty(spec), spec = cpmg_spec(); end
nS = numel(spec.ufls.f_Hz);

%% ------------------------------------------------------------------ %%
%  Initialise
%% ------------------------------------------------------------------ %%
if isempty(st)
    st = struct();
    st.tripped    = false(1,nS);
    st.armed      = false(1,nS);
    st.arm_time   = nan(1,nS);
    st.trip_time  = nan(1,nS);
    st.shed_total = 0;
    st.latching   = true;      % set false to allow automatic restoration
    st.min_df     = 0;
end

st.min_df = min(st.min_df, df);

%% ------------------------------------------------------------------ %%
%  Stage logic
%% ------------------------------------------------------------------ %%
df_reset = spec.hz2pu(spec.ufls.reset_Hz);

for k = 1:nS
    pickup = spec.ufls.df_pu(k);

    if ~st.tripped(k)
        if df <= pickup
            if ~st.armed(k)
                st.armed(k)    = true;
                st.arm_time(k) = t;
            elseif (t - st.arm_time(k)) >= spec.ufls.delay_s(k)
                st.tripped(k)   = true;
                st.trip_time(k) = t;
            end
        else
            % Frequency recovered before the operating time elapsed
            st.armed(k)    = false;
            st.arm_time(k) = NaN;
        end
    else
        if ~st.latching && df >= df_reset
            st.tripped(k)   = false;
            st.armed(k)     = false;
            st.arm_time(k)  = NaN;
        end
    end
end

%% ------------------------------------------------------------------ %%
%  Cumulative shed fraction
%% ------------------------------------------------------------------ %%
shed = sum(spec.ufls.shed_frac(st.tripped));
shed = min(max(shed,0),1);
st.shed_total = shed;

end

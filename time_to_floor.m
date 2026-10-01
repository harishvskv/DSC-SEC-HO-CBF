function T = time_to_floor(h, hdot, a_d, h_floor)
%TIME_TO_FLOOR  Predicted time for a barrier to fall from h to h_floor.
%
%   T = TIME_TO_FLOOR(h, hdot, a_d, h_floor) returns the first positive root of
%
%       h + hdot*t - 0.5*a_d*t^2 = h_floor
%
%   SINGLE SOURCE OF TRUTH.  This formula appears in compute_scatd (as the
%   SC-ATD bound), in estimate_contraction_rates (to extrapolate the fitted
%   decay) and in da_etm (as the predictive trigger).  It previously existed
%   as three private copies, one of which took the wrong root for negative
%   curvature and understated the tolerable duration.  Keeping one copy is the
%   only way that class of defect stays fixed.
%
%   VALID FOR BOTH SIGNS OF a_d
%     a_d > 0  concave, accelerating decline
%     a_d < 0  convex, DECELERATING decline, which is what this plant actually
%              exhibits: the governor ramps in against the GRC and arrests the
%              fall, so the fitted curvature is negative
%     a_d = 0  linear decline
%
%   In every case the first crossing is
%
%       T = [ hdot + sqrt( hdot^2 + 2*a_d*(h - h_floor) ) ] / a_d
%
%   For a_d < 0 with hdot < 0 this selects the smaller of the two positive
%   roots, which is the first crossing.  A negative discriminant means the
%   trajectory turns around before reaching the floor, so T = Inf.
%
%   RETURNS
%     0     when the barrier is already at or below the floor
%     Inf   when the floor is never reached
%
%   Vectorised over h and hdot when a_d and h_floor are scalars.
%
%   See also COMPUTE_SCATD, DA_ETM, ESTIMATE_CONTRACTION_RATES.

h    = h(:);
hdot = hdot(:);
if isscalar(h) && ~isscalar(hdot), h    = repmat(h,size(hdot)); end
if isscalar(hdot) && ~isscalar(h), hdot = repmat(hdot,size(h)); end

T = inf(size(h));

already = h <= h_floor;
T(already) = 0;

act = ~already;
if ~any(act), return; end

if abs(a_d) < eps
    lin = act & (hdot < 0);
    T(lin) = (h_floor - h(lin))./hdot(lin);
    return
end

disc = hdot.^2 + 2*a_d*(h - h_floor);
ok   = act & (disc >= 0);

Troot = (hdot(ok) + sqrt(disc(ok)))/a_d;
Troot(~isreal(Troot) | Troot < 0) = Inf;
T(ok) = Troot;
end

function [s, Ceff, info] = trapSamplerStep(s, pose, t)
%TRAPSAMPLERSTEP One ISR of the capacitance acquisition model.
%
%   [s, Ceff, info] = trapSamplerStep(s, pose, t)
%
% Call once per controller interrupt. Records the true capacitance at the
% current pose, then returns the four channels as the controller would
% actually see them: each delayed by its own acquisition age.
%
% INPUTS
%   s     sampler state from trapSamplerInit()
%   pose  6-by-1 [rx;ry;rz;thx;thy;thz], metres and RADIANS, at time t
%   t     [s] current time
%
% OUTPUTS
%   s     updated state
%   Ceff  [F] 4-by-1 DELAYED capacitances, as delivered to the controller.
%         NaN on every channel until the buffer has filled past maxAge.
%   info  .ageUsed  [s] 4-by-1 staleness of each returned value
%         .tSample  [s] 4-by-1 time each value actually refers to
%         .CeffNow  [F] 4-by-1 the true value at t, for comparison
%         .valid    logical, false while the buffer is still filling
%         .lagError [F] 4-by-1 CeffNow - Ceff, the error delay introduces
%
% The delayed values come from linear interpolation in the ring buffer, so
% the acquisition instants need not land on ISR boundaries. That matters
% because dtChannel (150 us) is not generally a multiple of Tisr.
%
% See also TRAPSAMPLERINIT, TRAPSENSING, TRAPCAPACITANCE.

% ---- true capacitance now ----
C    = trapCapacitance(s.plant, pose);
Cnow = trapSensing(s.plant, C);

% ---- push into the ring buffer ----
s.head = s.head + 1;
idx = mod(s.head-1, s.nBuf) + 1;
s.tBuf(idx)   = t;
s.cBuf(:,idx) = Cnow;

% ---- read each channel back at its own retarded time ----
Ceff    = nan(s.nCh,1);
tSample = t - s.age;
tOldest = min(s.tBuf(~isnan(s.tBuf)));
s.ready = ~isnan(tOldest) && (tOldest <= t - s.maxAge);

% interp1 needs at least two distinct abscissae. With a very short delay
% budget (maxAge below one ISR) the buffer can satisfy the age test while
% still holding a single sample, so guard on the count as well -- otherwise
% this throws MATLAB:griddedInterpolant:DegenerateGridErrId instead of
% simply reporting that the buffer is not ready yet.
have = ~isnan(s.tBuf);
s.ready = s.ready && (nnz(have) >= 2);

if s.ready
    tv = s.tBuf(have);
    [tv, ord] = sort(tv);
    for k = 1:s.nCh
        cv = s.cBuf(k, have);
        Ceff(k) = interp1(tv, cv(ord), tSample(k), 'linear', 'extrap');
    end
end

if nargout > 2
    info.ageUsed  = s.age;
    info.tSample  = tSample;
    info.CeffNow  = Cnow;
    info.valid    = s.ready;
    info.lagError = Cnow - Ceff;
end
end

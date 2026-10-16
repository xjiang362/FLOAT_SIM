function s = trapSamplerInit(plant, Tisr, varargin)
%TRAPSAMPLERINIT Set up the discrete capacitance acquisition model.
%
%   s = trapSamplerInit(plant)                       % all timing from plant.readout
%   s = trapSamplerInit(plant, Tisr)                 % override the ISR period
%   s = trapSamplerInit(plant, Tisr, 'dtChannel',80e-6, 'tI2C',200e-6)
%
% TUNING THE DELAY. The defaults live in plant.readout (see
% trapDeviceGeometry.m) so a sweep can change them in one place:
%
%   plant.readout.tI2C = 200e-6;        % try a faster bus
%   s = trapSamplerInit(plant);
%
% or override per call without touching the plant:
%
%   s = trapSamplerInit(plant, Tisr, 'tI2C', 200e-6);
%
% Models how the four quadrant readout pairs actually reach the
% controller. They are NOT sampled simultaneously: the front end converts
% them one at a time, then the whole set is shipped over I2C. So every
% channel arrives stale, and by a DIFFERENT amount.
%
%   t0           channel 1 (C12,   terminals 1,2)   converted
%   t0 +   dt    channel 2 (C56,   terminals 5,6)
%   t0 + 2*dt    channel 3 (C910,  terminals 9,10)
%   t0 + 3*dt    channel 4 (C1314, terminals 13,14)
%   t0 + 3*dt + tI2C                                 set delivered
%
% Age at delivery, channel k of n:
%       age(k) = tI2C + dtChannel*(n - k + 0.5)
% so with the defaults channel 1 is 1050 us old and channel 4 is 600 us
% old when the controller uses them. The SKEW between channels (450 us
% here) matters as much as the absolute lag: a pose estimator that treats
% the four as simultaneous will infer a spurious tilt whenever the plate
% is moving.
%
% INPUTS
%   plant   struct from trapDeviceGeometry()
%   Tisr    [s] controller interrupt period; the sampler is stepped once
%           per ISR and buffers at that resolution
%
% NAME/VALUE  (each defaults to the corresponding plant.readout field)
%   'dtChannel' [s] inter-channel conversion spacing
%   'tI2C'      [s] transfer delay after the last conversion
%
% OUTPUT
%   s   sampler state, to be passed to trapSamplerStep
%
% See also TRAPSAMPLERSTEP, TRAPSENSING.

if ~isfield(plant,'readout')
    error('trapSamplerInit:noReadout', ...
        'plant.readout is missing; regenerate the plant with trapDeviceGeometry.');
end
if nargin < 2 || isempty(Tisr), Tisr = plant.readout.Tisr; end

ip = inputParser;
ip.addParameter('dtChannel', plant.readout.dtChannel, @(v) isscalar(v) && v >= 0);
ip.addParameter('tI2C',      plant.readout.tI2C,      @(v) isscalar(v) && v >= 0);
ip.parse(varargin{:});

s.plant     = plant;
s.Tisr      = Tisr;
s.dtChannel = ip.Results.dtChannel;
s.tI2C      = ip.Results.tI2C;

s.nCh = size(plant.idx.sensePair,1);
k     = (1:s.nCh)';
s.age = s.tI2C + s.dtChannel*(s.nCh - k + 0.5);   % [s] per-channel staleness
s.maxAge = max(s.age);
s.skew   = max(s.age) - min(s.age);

% Ring buffer of past Ceff, at ISR resolution, long enough to reach back
% past the oldest channel with room for interpolation.
s.nBuf = ceil(s.maxAge/Tisr) + 3;
s.tBuf = nan(1, s.nBuf);
s.cBuf = nan(s.nCh, s.nBuf);
s.head = 0;          % number of samples written
s.ready = false;     % true once the buffer reaches back far enough
end


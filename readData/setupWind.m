%% setup wind data for forward model
%
% 

function [wind, avewind] = setupWind( mfname, dt, nsmooth, umin, nave)

% setupWind: reads wind.mat data which is prepared by readWind function and
% saved separately. Also sets cleans up the data, removes zero wind
% velocity and uses a laplacian operator to smooth the wind velocity. Also
% performs similar operations for wind direction.
%
% input :-
%   
% mfname - name of .mat file with the raw data.
% dt     - timespacing between data points.
% nsmooth- number of smoothing steps.
% umin   - minimum wind velocity below whihc wind is not detected.
%

% Default values:
if nargin < 2, dt = 600.0;  end  % time between wind data (600 s = 10 min).
if nargin < 3, nsmooth = 0; end  % number of smoothing steps. 
if nargin < 4, umin = 0.0;  end  % minimum wind velocity.
if nargin < 5, nave = 1;   end   % number of time steps to average 

UZero    = umin;        % zero wind value
kmph2mps = 1000/3600;  % converts km/h to m/s
deg2rad = 2*pi/360;
windshift= pi;         % wind blows FROM this direction

% load .mat file holding the data
load(mfname)

% setup time
wind.dt    = dt;   
wind.time  = [0:length(wind.vel)]' * dt;       % in seconds

% convert units
wind.vel = wind.vel* kmph2mps;
wind.dir = wind.dir* deg2rad;

%% filter wind
%
% set velocity below umin to umin. 
%

lowwind_indx = find( wind.vel < umin );
wind.vel(lowwind_indx) = umin;

wind.dir = mod(wind.dir+windshift, 2*pi); 

%% smoothing steps
%
% use an averaging operator with homogeneous Neumann boundary condition to
% smooth the data and preserve the mean. We perform the same operation on 
% sinh(wind.dir) to get rid of singularity in wind direction.
%

% create an averaging operator 
e = ones(length(wind.vel), 1);
Av = spdiags( [1/4*e, 1/2*e, 1/4*e], [-1, 0, 1], length(e), length(e) );
Av(1,1) = 3/4;
Av(end,end) = 3/4;

% perform smoothing 
for i = 1:nsmooth
   wind.vel = Av*wind.vel;
   wind.dir = asinh( Av*sinh(wind.dir) );
end

% average the wind to reduce data size
aveindx = [1:nave:length(wind.vel)];
avewind.vel = [];%zeros(1,length(aveindx));
avewind.dir = [];%zeros(1, length(aveindx));
avewind.time = [];% zeros(1, length(aveindx));
kk =1;

for k = aveindx(1:end-1)
    disp(k)
    disp(wind.vel(k:k+nave-1))
    avewind.vel(kk)  = (1/nave)*sum(wind.vel(k:k+nave-1));
    avewind.dir(kk)  = (1/nave)*sum(wind.dir(k:k+nave-1));
    avewind.time(kk) = wind.time(k);
    kk = kk+1;
end
avewind.n = length(avewind.vel);
wind.n = length(wind.vel);
end
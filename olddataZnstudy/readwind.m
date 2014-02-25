function wind = readwind( fname, dt, nsmooth, umin )
% READWIND: Read the wind data file with columns [date/time, raw speed,
%   speed, angle] and then filter the data, removing any NaNs or
%   other invalid values.  
%       
%   Input parameters:
%
%         fname - Excel file name
%            dt - time spacing between each data point in seconds,
%                 assumed constant (default: 600 s)
%       nsmooth - number of smoothing iterations (default: 0)
%       umin    - minimum wind velocity (default: 0)
%
%   Output parameter: a structure named 'wind' with fields
%
%         wind.time - time, starting at 0 (s)
%         wind.vel  - wind velocity (m/s)
%         wind.dir  - wind direction measured ccw from the North (rad)

% Default values:
if nargin < 2, dt = 600.0;  end  % time between wind data (600 s = 10 min).
if nargin < 3, nsmooth = 0; end  % number of smoothing steps. 
if nargin < 4, umin = 0.0;  end  % minimum wind velocity.

UZero    = 0.1;        % zero wind value
kmph2mps = 1000/3600;  % converts km/h to m/s
windshift= pi;         % wind blows FROM this direction

% Read the Excel file.  The data had to be rewritten to "Excel 4.0
% Worksheet" format (Matlab doesn't handle Excel files very well). 
warning( 'OFF', 'MATLAB:xlsread:Mode' );
[num,txt,raw] = xlsread( fname );
warning( 'ON',  'MATLAB:xlsread:Mode' );
len = size(num,1);
if size(raw,2) == 4, wformat = 1;       % format for original wind data
else                 wformat = 2; end   % format for "PartII" data

% The date column can be ignored because of the data/time format.
% So we assume all data is sampled at intervals of length 'dt'.
wind.dt    = dt;   
wind.time  = [0:len-1]' * dt;       % in seconds

if wformat == 1, 
  wind.vel   = num(:,3) * kmph2mps;   % convert km/h to m/s
  wind.dir   = num(:,4) * (2*pi/360); % in radians
else
  wind.vel   = num(:,2) * kmph2mps;   % convert km/h to m/s
  wind.dir   = num(:,7) * (2*pi/360); % in radians
end

% Make sure wind lies in the interval [0,2*pi].  Also,
% shift the wind direction so the wind blows TO the
% given direction and not FROM (as it is in the original 
% data file.
wind.dir = mod(wind.dir+windshift, 2*pi);           
  
%----------------------------------------------------------------------
% Filter the velocity:
%   1) For any wind velocities less than umin, set them equal to umin.
%      At the same time, the wind directions for these values are set
%      to NaN so that they are interpolated from neighbouring values
%      (with the reasoning that anemometers in calm winds give an
%      unreliable measurement of direction).
%   2) Replace all NaN's appearing in the wind direction with
%      values interpolated linearly from neighbouring non-NaN
%      values. 
%   3) Introduce a Laplacian smoothing into the direction.
%      ??? ALSO SMOOTH THE SPEED ???
%----------------------------------------------------------------------

% 1) Set the minimum wind velocity to umin.  Empty values are read in
%    as NaN by XLSREAD.
ii  = find(isnan(wind.vel));
wind.vel(ii) = UZero;
ii2 = find(wind.vel < umin);
wind.vel(ii2) = UZero;
wind.dir(ii2) = NaN;
fracmin = (length(ii) + length(ii2)) / length(wind.vel);
fprintf( 1, 'Fraction of wind measurements less than umin=%f is %f\n', umin, fracmin );

% 2) Use linear interpolation to fill in any NaN values in 
%    the wind direction.
if isnan(wind.dir(1)),
  wind.dir(1) = wind.dir(find(~isnan(wind.dir), 1));
end
if isnan(wind.dir(end)), 
  wind.dir(end) = wind.dir(find(~isnan(wind.dir), 1, 'last'));
end
  
i = 2;
while i < len,
  j = find(isnan(wind.dir(i:end)),1)+i-1;
  jj= find(~isnan(wind.dir(j:end)),1)+j-1;
  if ~isempty(jj),
    ddir = (wind.dir(jj)-wind.dir(j-1)) / (jj-j+1);
    wind.dir(j:jj-1) = [1:(jj-j)]' * ddir + wind.dir(j-1);
  end
  i = jj;
end
  
% 3) Apply a simple "averaging" smoother to the direction values. 
%    Wind speed is NOT smoothed.
%
%    NOTE: This still doesn't take care of wrap-around from 0 to 360
%          degrees. Should we do something like average the
%          sine/cosine/tangent of the angle instead? 
for i = 1 : nsmooth,
  wind.dir(2:end-1) = (wind.dir(1:end-2) + 2*wind.dir(2:end-1) + wind.dir(3:end))/4;
  wind.vel(2:end-1) = (wind.vel(1:end-2) + 2*wind.vel(2:end-1) + wind.vel(3:end))/4;
end

% END.

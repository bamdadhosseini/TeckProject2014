% FORWARD: Solve the forward atmospheric dispersion problem.  That
%    is, given the source emission rates (in mol/s), calculate the
%    amount of Zn deposited in each receptor (in kg).

addpath('.././readData/');
addpath('.././data/');

tic
clear all

setparams;

% Add noise to the data (if necessary).
%dnoise = 0.0;
%source.Q = source.Q .* (1 + dnoise*(1-2*rand([1,4])));  % add noise

% setup solver parameters
setSolver;
source.x = source.x';
source.y = source.y';
source.z = source.z';
recept.x = recept.x';
recept.y = recept.y';
recept.z = recept.z';
source.Q = ones(size(source.Q));

% wind data and time are stored in a separate .mat file to accelerate
% loading. These are extracted from original excel files under './data'
wind = setupWind( 'windData.mat', 600, 5, umin);

if isanim, figure(2), clf, end
fprintf( 1, 'Accumulating deposition data (tskip=%d) ...\n', tskip );
dt  = (wind.time(2) - wind.time(1)) * tskip;
dep = zeros(source.n, recept.n);
dep2 = zeros(size(xmesh,1), size(xmesh,2), source.n);
nwind = length(wind.dir);

for k = [1 : tskip : nwind],

  k2 = min( k+tskip-1, nwind );
  if IsAvgWind,
    thetaavg = mean( wind.dir(k:k2) );
    Vwavg    = mean( wind.vel(k:k2) );
  else
    thetaavg = wind.dir(k);
    Vwavg    = wind.vel(k);
  end
  
  % in case we want animations
  animate;

  if floor((k-tskip)/500) ~= floor(k/500), fprintf( 1, ' %5d ..', k ), end
  theta = +(pi/2 - thetaavg);    % ** CHECK THIS! **
  Vw    = Vwavg;
  warning( 'OFF', 'MATLAB:divideByZero' );

%     %% First, accumulate depositions at the receptor sites.
%     xxx = (recept.x - source.x(i))*cos(theta) + (recept.y - source.y(i))*sin(theta);
%     yyy =-(recept.x - source.x(i))*sin(theta) + (recept.y - source.y(i))*cos(theta);
%     C = ermak( xxx, yyy, recept.z, source.z(i), source.Q(i), ...
%                Vw, Vspbo, Vdpbo, stabclass );
%     % C is a concentration in mol/m^3, convert to deposition in kg ...
%     dep = dep + (A * dt * Mpb * Vdpbo) * C; 

  dep = getDepAtRecept( dep,recept, source, theta, Vw, Vspbo, Vdpbo, stabclass ...
                     ,A, Mpb, dt);
  %dep = sum(dep,1);
  if isplot,
  % Next, accumulate depositions in all cells on the mesh.
  dep2 = getDepDomain( dep2, xmesh, ymesh, zmesh, dx, dy, source, theta,...
                            Vw, Vspbo,Vdpbo, stabclass, Mpb, dt); 
  end
  warning( 'ON', 'MATLAB:divideByZero' );
  
end

% sum over deposition from all sources
dep = sum(dep,1);
dep2 = sum(dep2, 3);

fprintf( 1, '\nDeposition in each receptor (mg):\n' );
disp(dep(1:nrec)*1e6)
fprintf( 1, 'Total deposited in all receptors (mg):\n' );
disp(sum(dep(1:nrec))*1e6)

receptDepPlot;

if isplot,
   domainDepPlot;
end

save 'forsave.mat' 
toc

rmpath('.././readData/');
rmpath('.././data/');
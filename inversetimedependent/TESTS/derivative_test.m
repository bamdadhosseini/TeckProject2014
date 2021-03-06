%% Derivative test 
%
% After the adjoint test we need to perform the derivative test. The
% functional consists of three terms representing the misfit and a fourth
% regularization term. To keep things simple we choose the regularization
% to be the laplacian operator. So the functional is 
%
% J = 1/2* ( \sum_i \gamma_i || M_iEQ - data_i ||_2^2 ) 
%            + 1/2 * \sum_j|| LQ_j ||_2^2
%
% where \gamma_i is the weight for each type of data (xact, TSP, dustfall,
% ...), M_i is the corresponding measurement operator and E is the forward 
% map for concentrations. L is the regularization term.
%
%

clear all
close all
clc

%profile on

addpath('.././readData/');
addpath('.././data/inversedata');
addpath('.././problemsetup');
addpath('.././util');
addpath('.././forwardtimedependent');

setSolver;
[realwind, wind] = setupWind('windData.mat', 600, 5, umin, 6);
nwind = wind.n;
dt  = (wind.time(2) - wind.time(1)) * tskip; % assuming wind is measured on 
                                             % constant intervals
disp(dt)
setparams;                         
                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% construct a random matrix for source emissions (no real need for a smooth
% one but it won't hurt)

source.Q = zeros(source.n, nwind);
sourcevariance = 2e-5;
for i =1:source.n
   source.Q(i,:) = getSampleSmoothSource(nwind, sourcevariance)*1e-1; % mol/s 
end
figure(1);
plot(wind.time(1:nwind), source.Q);
xlabel('time (matlab numeric)');
ylabel('mol/s');

% construct a second set of random functions as perturbation to the
% sources to use in the derivative test.

deltaQ = zeros(source.n, nwind);
pertvariance = max(max(source.Q));
e = ones(nwind,1);
L = spdiags([-e, 2*e, -e], [-1,0,1], nwind, nwind);
for i =1:source.n
   deltaQ(i,:) = L\(pertvariance*(randn(nwind,1)))*1e-4; % mol/s 
end
figure(2);
plot(wind.time(1:nwind), deltaQ);
xlabel('time (matlab numeric)');
ylabel('mol/s');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------------------');
disp('constructing synthetic data');

noiselevel = 0.3;
% solve the forward problem once to get synthetic data 
tic
C = getConcentrationAtReceptors(sensor, source, wind, Vspbo, Vdpbo, stabclass, tskip);
sensor = getAccumulatedConcentration( C, sensor, 1:length(sensor));
toc
% take measurements and construct data structure
% measurement.dustfall = getDepositionAtDustfallJar ( sensor(sensorindex.dustfall), A, dt, Mpb, Vdpbo );
% 
% data.dustfall = num2cell(measurement.dustfall + ...
%     0.3*max(max(measurement.dustfall))*rand(size(measurement.dustfall)));
% [sensor(sensorindex.dustfall).data]  = deal(data.dustfall{:});
% 
% measurement.xact     =  Mpb*getAverageConcentration( sensor(sensorindex.xact) )';
% data.xact = num2cell(measurement.xact + ...
%     0.3*max(max(measurement.xact))*rand(size(measurement.xact)), 1);
% [sensor(sensorindex.xact).data]     = deal(data.xact{:});
% 
% 
% measurement.TSP      =  Mpb*getAverageConcentration( sensor(sensorindex.TSP));
% data.TSP = num2cell(measurement.TSP + ...
%     0.3*max(max(measurement.TSP))*rand(size(measurement.TSP)), 2);
% [sensor(sensorindex.TSP).data]     = deal(data.TSP{:});
% 
% measurement.PM10     =  Mpb*getAverageConcentration( sensor(sensorindex.PM10));
% data.PM10 = num2cell(measurement.PM10 + ...
%     0.3*max(max(measurement.PM10))*rand(size(measurement.PM10)),2);
% [sensor(sensorindex.PM10).data]     = deal(data.PM10{:});

sensor = getMeasurementAtSensor( sensor );
data   = arrayfun( @(s) abs(s.measurement + ...
        noiselevel*(mean(s.measurement))*(rand(size(s.measurement))-1/2)),...
        sensor, 'UniformOutput', 0);
% distribute data to sensors 
[sensor(:).data] = deal( data{:} );
clear data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% derviative test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('');
disp('----------------------------');
disp('Starting Derivative Test');
disp('');
%
% how much do we trust each sensor?
%
% just random here
%

trustweight(sensorindex.dustfall) = 1/3;
trustweight(sensorindex.xact) = 1/6;
trustweight(sensorindex.TSP) = 1/6;
trustweight(sensorindex.PM10) = 1/3;
trustweight = num2cell(trustweight);


[sensor(:).weight] = deal(trustweight{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compute misfits
%
sensor = getMeasurementAtSensor( sensor );
sensor = getMisfitAtSensor( sensor );

% 
% compute the Tikhonov functional at this point
%
alpha = 1;
R = defineRegularization( nwind, dt );

J = evalTikhonovFunctional( source, sensor, R, alpha);
disp(['J = ', num2str(J)]);
%
% compute the derivative using adjoints
%
tic
adjsource = getAdjointSourceForSensor( sensor, size(C) );
adjsol = solveAdjointConcentrationAtReceptors( sensor,...
    source, wind, Vspbo, Vdpbo, stabclass, ...
    tskip, adjsource);
%
% derivative of regularization
%
gradR = source.Q*(R'*R);

gradJ = adjsol(:) + alpha*gradR(:);
toc
%
% perform an adjoint test to make sure
%
disp('adjoint test');
%
% stack all measurements on top of each other
%
allmeas = stackOnTop( sensor, 'measurement' );
%
% the WW vector carries the effect of trust wirghts, this is only 
% used for the adjoint test.
%
WW = arrayfun( @(s) ones(numel(s.measurement),1).*s.weight,...
    sensor, 'UniformOutput', 0);
WW = cat(1, WW{:});
allmisfit = stackOnTop( sensor, 'misfit' );
disp( ['result = ', num2str(source.Q(:)'*adjsol(:) -  (WW.*allmeas)'*allmisfit)] );
%pause

eps = 5.^(0:-1:-10);
cntr = 1;

disp('');
disp('----------------------------------');
disp(' iterating over perturbation size ');
disp('----------------------------------');
disp('');
pertsource = source;
for pertSize = eps
   
    % perturb the source term
    pertsource.Q = source.Q + pertSize*deltaQ;
    
    % solve forward problem and compute data misfit for perturbed data 
    tic
    C = getConcentrationAtReceptors(sensor, pertsource, wind, Vspbo, Vdpbo, stabclass, tskip);
    sensor = getAccumulatedConcentration( C, sensor, 1:length(sensor));
    sensor = getMeasurementAtSensor( sensor );
    sensor = getMisfitAtSensor( sensor );
    
    % compute the Tikhonov functional at new point
    Jpert = evalTikhonovFunctional( pertsource, sensor, R, alpha);
    toc
    disp(['Jpert = ', num2str(Jpert)]);
    % perform derivative test
    error1(cntr) = abs( Jpert - J ); 
    error2(cntr) = abs( Jpert - J - pertSize*(gradJ'*deltaQ(:)) );
    
    disp('');
    disp([num2str(cntr), ' eps = ', num2str(pertSize),...
        '           ', num2str(error1(cntr)), ', ', num2str(error2(cntr))]);
    cntr = cntr+1;
end

figure(3);
loglog( eps, error1, 'r', eps, error2, 'b' );
xlabel( 'perturbation size \epsilon' );
ylabel( 'error' );

legend('| J(Q + \delta Q) - J(Q) |', ...
       '| J(Q + \delta Q) - J(Q) - \epsilon \nabla J \delta Q |'); 

toc

rmpath('.././readData/');
rmpath('.././data/inversedata');
rmpath('.././problemsetup');
rmpath('.././util');
rmpath('.././forwardtimedependent');
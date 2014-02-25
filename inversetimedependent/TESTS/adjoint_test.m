%% adjoint tests for timedep inverse problem

% the inverse problem consists of minimizing three functionals and we 
% have to identify the adjoint of each term. This script performs the 
% test for each one of these terms as well as identifying the measurement
% operators and their adjoints. This is not trivial since we cannot store
% the whole forwardmap or the entire measurement operator. Instead we 
% apply it explicitly. For example, the forward map for accumulating 
% deposition of all sources is a row of block identitiy matrices which is 
% very large for each time step. It's adjoint is a column vector of
% identity matricies again. We dont have enough memory to store them so 
% we identify what they do and simply apply that in each time step. 
%


clear all
close all
clc

profile on

% setup the problem first


addpath('.././readData/');
%addpath('.././data/');
addpath('.././data/inversedata');
addpath('.././problemsetup');
addpath('.././util');
addpath('.././forwardtimedependent');

setSolver;
wind = setupWind('windData.mat', 600, 5, umin);
nwind = wind.n;
dt  = (wind.time(2) - wind.time(1)) * tskip; % assuming wind is measured on 
                                             % constant intervals
setparams;                         
                                             
% construct a random matrix for source emissions (no real need for a smooth
% one but it won't hurt)

source.Q = zeros(source.n, nwind);
sourcevariance = 2e-5;
for i =1:source.n
   source.Q(i,:) = getSampleSmoothSource(nwind, sourcevariance)*1e-3; % mol/s 
end
figure(1);
plot(wind.time(1:end-1), source.Q);
xlabel('time (matlab numeric)');
ylabel('mol/s');

%% Foward map 
% Apply the forward map at each time step to get concentration and then
% calculate depositions
tic 
disp('Solving Forward problem');
C = getConcentrationAtReceptors(sensor, source, wind, Vspbo, Vdpbo, stabclass, tskip);

%% accumulated concentration
% 
% all sensors work with accumulated concentration of solution at 
% various time scales. This functions computes these values for all
% kinds of sensors.

sensor = getAccumulatedConcentration( C, sensor, 1:length(sensor));
toc
%% Solve forward
%
% apply measurement operator M_d (which is simply computing depositions and
% averages) we normally wouldn't store these separately but this helps
% clear things up for the adjoint test.
%
tic
disp('Taking measurements');

measurement.dustfall =  getDepositionAtDustfallJar ( sensor(sensorindex.dustfall), A, dt, Mpb, Vdpbo );

measurement.xact     =  Mpb*getAverageConcentration( sensor(sensorindex.xact) );

measurement.TSP      =  Mpb*getAverageConcentration( sensor(sensorindex.TSP));

measurement.PM10     =  Mpb*getAverageConcentration( sensor(sensorindex.PM10));

measurement.all      =  [ measurement.dustfall; ...
                          measurement.xact'; ...
                          reshape(measurement.TSP', numel(measurement.TSP), 1);
                          reshape(measurement.PM10', numel(measurement.PM10), 1)
                          ];
toc


%
% Adjoint problem
%

%
% construct random misfit at each sensor
%

for iSens = 1:recept.n
    sensor(iSens).misfit = rand(size(sensor(iSens).accum_concentration));
end

%
% construct adjoint source terms based on the misfit data
% again we dont really need them as separate sources but this 
% gives more insight into the solution. We can test things separately
% now.
%
tic
disp('Solving the Adjoint problem');
adjsource = getAdjointSourceForSensor(sensor, size(C));

% scale adjoint sources properly 
% this is not required in real problems since misfit is already 
% scaled properly.

adjsource(sensorindex.dustfall,:) = adjsource(sensorindex.dustfall,:).*...
    (A * dt * Mpb * Vdpbo);

adjsource(sensorindex.xact,:) = adjsource(sensorindex.xact,:).*...
    (Mpb*1/sensor(31).accum_period);

adjsource(sensorindex.TSP,:) = adjsource(sensorindex.TSP,:).*...
    (Mpb*1/sensor(32).accum_period);


adjsource(sensorindex.PM10,:) = adjsource(sensorindex.PM10,:).*...
    (Mpb*1/sensor(34).accum_period);

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

toc

% adjoint test for different data

% deposition
tic
disp('')
disp('---------------------');
disp('Adjoint test for dustfalljar');
adjsolution = solveAdjointConcentrationAtReceptors( sensor(sensorindex.dustfall), source, wind, Vspbo, Vdpbo, stabclass, ...
    tskip, adjsource(sensorindex.dustfall,:));

% adjoint test 
disp( source.Q(:)'*adjsolution(:) - ...
    measurement.dustfall(:)'*[sensor(sensorindex.dustfall).misfit]');
toc

% Xact
tic
disp('')
disp('---------------------');
disp('Adjoint test for Xact');
adjsolution = solveAdjointConcentrationAtReceptors( sensor(sensorindex.xact),...
    source, wind, Vspbo, Vdpbo, stabclass, ...
    tskip, adjsource(sensorindex.xact,:));

% adjoint test 
disp( source.Q(:)'*adjsolution(:) - ...
    measurement.xact(:)'*[sensor(sensorindex.xact).misfit]);
toc

% TSP
tic
disp('')
disp('---------------------');
disp('Adjoint test for TSP');
adjsolution = solveAdjointConcentrationAtReceptors( sensor(sensorindex.TSP),...
    source, wind, Vspbo, Vdpbo, stabclass, ...
    tskip, adjsource(sensorindex.TSP,:));

% adjoint test 
TT = measurement.TSP';

disp( source.Q(:)'*adjsolution(:) - ...
    TT(:)'*[sensor(sensorindex.TSP).misfit]');
toc

% PM10
tic
disp('')
disp('---------------------');
disp('Adjoint test for PM10');
adjsolution = solveAdjointConcentrationAtReceptors( sensor(sensorindex.PM10),...
    source, wind, Vspbo, Vdpbo, stabclass, ...
    tskip, adjsource(sensorindex.PM10,:));

% adjoint test 
PP = measurement.PM10';
disp( source.Q(:)'*adjsolution(:) - ...
    PP(:)'*[sensor(sensorindex.PM10).misfit]');
toc

% All at once
tic
disp('')
disp('---------------------');
disp('Adjoint test all at once');
adjsolution = solveAdjointConcentrationAtReceptors( sensor,...
    source, wind, Vspbo, Vdpbo, stabclass, ...
    tskip, adjsource);

% adjoint test
sensor(31).misfit = sensor(31).misfit'; % im too lazy to fix this right now
disp( source.Q(:)'*adjsolution(:) - ...
    measurement.all'*[sensor(:).misfit]');
toc

profile viewer

rmpath('.././readData/');
rmpath('.././data/inversedata');
rmpath('.././problemsetup');
rmpath('.././util');
rmpath('.././forwardtimedependent');


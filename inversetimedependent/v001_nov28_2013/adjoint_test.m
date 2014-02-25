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

% setup the problem first


addpath('.././readData/');
addpath('.././data/');
addpath('.././data/inversedata');
addpath('.././forwardsolve');
addpath('.././util');
addpath('.././forwardtimedependent');

setparams;
setSolver;
wind = setupWind('windData.mat', 600, 5, umin);
nwind = wind.n;
dt  = (wind.time(2) - wind.time(1)) * tskip; % assuming wind is measured on 
                                             % constant intervals
                         
                                             
%dxact = readXact('xact_data.xls');
%dTSP = readTSP('TSP_data.xls');

%%% Adjoint test for deposition %%%%%%%%%%%%%%%
% given a random vector of time dependent sources Q  and random vector of 
% deposition misfit d then if E is the forward map for concentrations and 
% M_d being the map that measures the depositions then for the adjoint test
% we want
%
% ( M_d * E * Q , d ) - ( Q, E' * M_d' * d ) = 0
%
% M_d * E * Q  = deposition at each receptor.
% M_d' * d     = a tall vector of d being repeated at each time 
% E' at each time step is simply the transpose of the forward map at that 
% time. 

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

C = getConcentrationAtReceptors(recept, source, wind, Vspbo, Vdpbo, stabclass, tskip);

%% deposition
% apply measurement operator M_d (which is simply computing depositions)
M_dEQ = getDepositionAtReceptors( C, A, dt, Mpb, Vdpbo );

% Adjoint map
%

% construct a random vector of misfit ( this does not have a physical
% meaning it will only be used for adjoint test )
d = rand(recept.n, 1);

% now apply the M_d' to turn d into a vector of constant size for all time
% steps 
M_dTd = (A * dt * Mpb * Vdpbo)*repmat(d, 1, wind.n); % this serves as the source term in the adjoint

M_dTETd = solveAdjointConcentrationAtReceptors( recept, source, wind, Vspbo, Vdpbo, stabclass, ...
    tskip, M_dTd);

toc

% adjoint test 
disp(' Adjoint test for deposition measurement ');
disp(abs(M_dEQ'*d - source.Q(:)'*M_dTETd(:)));

%% Time dependent measurements
% time averages of concentrations ( for Xact, TSP and PM10 data )
%
% here the measurement operator starts at some given time step t_start 
% and until it reaches t_end it takes an average of the concentration 
% at receptor locations every m time steps.

TSP(1).start_indx = 101:300:3001;   % index of time step when measurement starts 
                                    % note that if the sensor starts at
                                    % indx = 100 for example, then it
                                    % measures the 101th otherwise
                                    % measurement period would not match.
TSP(1).end_indx   = 150:300:3050;
TSP(2).start_indx = 1;
TSP(2).end_indx   = 4608;
TSP(1).ave_period = 10;     % number of averaging time steps 
TSP(2).ave_period = floor(4609/2);

%disp(['test ', num2str(sum(C(2,:))/numel(C(2,:)))]);

% take measurements given the start and finish times 
% just to keep things simple we assume the TSPs are 
% placed at the first and second receptors.
% this completes the forward solve for the time dependent measurements
TSP = getAverageConcentration( C(1:2,:), TSP, Mpb, dt ); % say the last two rows are the TSPs

% now apply the adjoint of the measurement operator. This Computes the
% misfit and turns it into a source term which is constant during the
% measurement intervals for the adjoint equation. So it is a piecwise 
% constant source in time.

% construct a random vector of misfit for each sensor 
TSP(1).misfit = rand(size(TSP(1).forwardmeasurement));
TSP(2).misfit = rand(size(TSP(2).forwardmeasurement));

% compute the source term of the adjoint equation
M_STd = getAdjointSourceForSensor(TSP, size(C(1:2,:)));

% compute solution of the adjoint equation
%sensor.n = 2
%M_STESd = getAdjointSolutionForSensor( recept(1,:), source, wind, Vspbo, Vdpbo, stabclass, ...
%    tskip, M_STd )

rmpath('.././readData/');
rmpath('.././data/');
rmpath('.././data/inversedata');
rmpath('.././forwardsolve');
rmpath('.././util');
rmpath('.././forwardtimedependent');


%% Inverse Descent with constraints

clear all
close all
clc

%profile on

addpath('.././readData/');
addpath('.././data/inversedata');
addpath('.././problemsetup');
addpath('.././util');
addpath('.././forwardtimedependent');

global verbosity;
verbosity = 1;
scrsz = get(0,'ScreenSize');
h = figure('Position',[50 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])

%% setup parameters and read wind data 

setSolver;
[fullwind, wind] = setupWind('windData.mat', 600, 5, umin, 6);
nwind = wind.n;
dt  = (wind.time(2) - wind.time(1)) * tskip; % assuming wind is measured on 
                                             % constant intervals
setparams_john; 

%% get data
%
% so far it is synthetic to test the algorithm
noiselevel = 0.1; % 10 percent
originalsource = source; % we keep track of the data generator 
                         % for later comparison
[sensor, originalsource] = getSyntheticData( originalsource, sensor,...
                            wind, Vspbo, Vdpbo,...
                            stabclass, tskip, noiselevel );  

% we want the code to target the averages of the sources
%
targetQ = repmat( mean(originalsource.Q, 2), 1, size(originalsource.Q, 2));

%% setup regularization
%
alpha = 1e-2;
R = defineRegularization( wind.n, dt, 'periodic' );           

%% Descent algortihm
%

optmethod ='GPRN';
numDescentSteps = 100;

% start at zero (just a compromise)
source.Q = (1)*ones(source.n, wind.n);

% line search parameters
Armijomaxiter = 20;
tau = 1/2; % linesearch parameter
c1 = 1e-4; % linesearch parameter
gammak = 1;

% setup box constraints 
l = zeros(size(source.Q));
u = (1e3)*ones(size(source.Q));
%u(end/2:end) = 2;

% aux variables
mem = 4; % number of previous gradient and solutions to remember
s = zeros(source.n*wind.n, mem);
y = zeros(source.n*wind.n, mem);
gradJold = zeros(source.n, wind.n);
Qold     = zeros(source.n, wind.n);

%%%%%%%%%%%%%%%%%%%%%%
%%%%% END OF SETUP %%%
%%%%%%%%%%%%%%%%%%%%%%

for i = 1:numDescentSteps

    disp('')
    disp(['step: ', num2str(i)]);
    tic
    
    if (i==1)
        % compute the objective at initial guess 
        % to start the algorithm
        [J(i), C, sensor] = computeObjective(sensor, source, wind, ...
                    Vspbo, Vdpbo, tskip, stabclass, R, alpha, targetQ);
    end
    
    %% get descent direction
    %
    if (i>1)
        gradJold = gradJ; 
    end;
    
    gradJ = computeGradient(sensor, source, wind, ...
                    Vspbo, Vdpbo, tskip, stabclass, R, alpha, C, targetQ);
    %[descentDir] = getDescentDir( gradJ, optmethod, R, i, l, u);
    [descentDir] = getDescentDirConstraint( gradJ, optmethod, R, i, l, u, ...
    sensor, source, wind, Vspbo, Vdpbo,...
        tskip, stabclass, alpha, targetQ, J, c1, tau, gammak, Armijomaxiter);
    
    %% line search using armijo
    %
    [Qk, gammak, Jk, C, sensor] = armijoLineSearch(sensor, source, wind, Vspbo, Vdpbo,...
        tskip, stabclass, R, alpha, targetQ, descentDir, J, gradJ, c1, tau, gammak, Armijomaxiter,i, l, u );
    %% update
    %
    Qold = source.Q;
    source.Q = Qk;
    J(i+1) = Jk;
    
    
    % plot progress
    %
    figure(h);
    a1 = subplot(2,2,1); plot(wind.time(1:wind.n), originalsource.Q);
    xlabel('time'); ylabel('Q'); title('Target');
    
    a2 = subplot(2,2,2); plot(wind.time(1:wind.n), source.Q);
    xlabel('time'); ylabel('Q'); title(['Solution at step ', num2str(i)]);
    
    %linkaxes([a1 a2], 'y');
    %set(a1,'Ylimmode','auto'); 
    
    subplot(2,2,3, 'Yscale', 'log'), hold on 
    loglog(J, '-s');
    xlabel('step'); ylabel('J'); title('Objective Functional');
    
    subplot(2,2,4), plot(wind.time(1:wind.n), descentDir);
    xlabel('time'); ylabel('\nabla J'); title('Current Direction');
    

drawnow;
    
end
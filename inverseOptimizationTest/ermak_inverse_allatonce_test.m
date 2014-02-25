%% All at once inverse solve

%% ermak_inverse_test


clear all
%clc
close all

tic

setparams;
setSolver;

%% create synthetic data

synthsource = source;
synthsource.Qzn = [0.1 1 0.01 0.01];
[dep, dep2] = forwardErmak(synthsource, recept, wind, xmesh,...
    ymesh, zmesh, dx, dy, dz, dt, tskip, nwind, Vszns, Vdzns, Mzn, stabclass,A,...
    IsAvgWind, IsOverDomain);

totaldep = sum(dep,1); % sum over deposition from all sources to get
                       % the model output
eps = 1e-2;
synthdep = totaldep + eps*rand(1,recept.n)*max(totaldep);
d = synthdep;
 
%%
%% Solve Inverse %%%%%%%%%%%%%%
%%%  uses steepest descent %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up other solver parameters. 
nx = 10;
ny = nx;
nz = 10;
xlim = [-300,1400];
ylim = [-600,1400];
zlim = [0,100];
dx = (xlim(2)-xlim(1)) / nx;
dy = (ylim(2)-ylim(1)) / ny;
dz = (zlim(2)-zlim(1)) / nz;
x0 = [xlim(1):dx:xlim(2)]; % distance along the wind direction (m)
y0 = [ylim(1):dy:ylim(2)]; % cross-wind distance (m)
z0 = [zlim(1):dz:zlim(2)]; % vertical distance (m)
[xmesh,ymesh] = meshgrid(x0,y0);
zmesh = 0; 

setOptimizer;

%% Precomputations
%
% Compute deposition of each source with unit emission 
% rate and store this in an [S X M] matrix where S is 
% the number of sources and M is the number of receptors
% This Matrix is Called "E" for the Ermak solver.
%
unitSource = source;
unitSource.Qzn = unitSource.Qzn./unitSource.Qzn; % lazy way to create my unit sources

% Compute the matrix E
[E, E2] = forwardErmak(unitSource, recept, wind, xmesh,...
    ymesh, zmesh, dx, dy, dz, dt, tskip, nwind, Vszns, Vdzns, Mzn, stabclass,A,...
    IsAvgWind, IsOverDomain); 

%% Steepest Descent
%
% update solution as 
% 
% Q_{k+1} = Q_k + \beta_k P_k
%
% where \beta_k is the step size and 
% P_k is the search direction. In Steepest
% Descent we take P_k to be -eye(numel(gradJ))gradJ
% and also choose \beta_k = gradJ'*gradJ/(gradJ' Q gradJ)

alpha = 1e-8;
disp('All at once');
disp('');

%% Solve for Q

% solve for Q in one step by setting gradJ = 0
Q = (E*E' + alpha*eye(source.n))\(E*d');

disp('===============================================');
disp(['Original Source  :', num2str(synthsource.Qzn)]);
disp(['Solution         :', num2str(Q')]);
disp('===============================================');

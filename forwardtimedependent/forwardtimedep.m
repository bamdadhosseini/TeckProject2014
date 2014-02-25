%% script for solving the forward problem with variable sources 

clear all
close all
clc

% setup the problem first


addpath('.././readData/');
addpath('.././data/');
addpath('.././data/inversedata');
addpath('.././forwardsolve');
addpath('.././util');

tic

setparams;
setSolver;
wind = setupWind('windData.mat', 600, 5, umin);
nwind = wind.n;
dt  = (wind.time(2) - wind.time(1)) * tskip; % assuming wind is measured on 
                                             % constant intervals
                         
                                             
%dxact = readXact('xact_data.xls');
%dTSP = readTSP('TSP_data.xls');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct a sample time varying set of sources

source.Q = zeros(source.n, nwind);
sourcevariance = 2e-5;
for i =1:source.n
   source.Q(i,:) = getSampleSmoothSource(nwind, sourcevariance)*1e-3; % mol/s 
end
figure(1);
plot(wind.time(1:end-1), source.Q);
xlabel('time (matlab numeric)');
ylabel('mol/s');

% if you are testing with previous dep code
%source.Q = ones(size(source.Q));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% construct forward map for concentration 
%
% Matrix E such that 
% 
%  E*vec(Q) = vec(C)
% 
% gives concentration of the pollutant at each receptor.
%

% correct wind dir 
wind.dir= +(pi/2 - wind.dir);
C = getConcentrationAtReceptors( recept, source, wind, Vspbo, Vdpbo, stabclass, tskip);
figure(2)
plot(wind.time(1:end-1), C*Mpb*1e6);
xlabel('time');
ylabel('local concentration mg/m^3');

% compute deposition at receptors
d = (A * dt * Mpb * Vdpbo)*sum(C, 2);

toc 

disp( 'deposition at receptors')
disp( d' );
figure(3)
bar( d*1e6 );
xlabel( 'receptor');
ylabel('deposition in mg');
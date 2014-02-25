%% Inverse solution using Tikhonov regularization
%
% solves the normal equation :
%
% Q(EE' + alpha*I)  = data*E'
%
% where E is the forward solution map for the ermak solution and 
% Q is the row vector of  emissions rate, alpha is the regularization
% parameter and  data is measured deposition. The solution is the
% minimizer of 
%
% J_alpha = || QE - data||_2^2 + alpha*||Q||_2^2
%
%

clear all
close all
clc

addpath('.././readData/');
addpath('.././data/');
addpath('.././data/inversedata');
addpath('.././forwardsolve');
addpath('.././util');

tic

setparams;
setSolver;
wind = setupWind('windData.mat', 600, 5, umin);
nwind = length(wind.vel);
dt  = (wind.time(2) - wind.time(1)) * tskip; % assuming wind is measured on 
                                             % constant intervals
                                        
% QQ = 10*rand(source.n,1);
% source.Q = QQ;
% [E, E2] = forwardErmak( source, recept, wind, ...
%                         xmesh, ymesh, zmesh, dx, dy, dz, dt, tskip, nwind,...
%                         Vspbo, Vdpbo, Mpb, stabclass, A, ...
%                         IsAvgWind, IsOverDomain );
% noiselev = 1e-3;
% d = sum(E, 1) + noiselev*rand(1,size(E,2));

%% read lead dep data

d = readdustfall('data_Pb.xls');

% reset receptors, leave out those that dont have available data 
recept = clearReceptors( recept, d );

d = d(find(d ~=0))'*1e-6; % convert mg to kg of deposition

% add some random noise?
%
% make sure the noise does not result in negative data
noiselevel= 0.3;
while true
   
   variation = (rand(size(d))-1/2)*noiselevel*mean(d);
   if (d + variation > 0)
       d = d + variation;
       break;
   end
end
%d = d + (rand(size(d))-1/2)*noiselevel*mean(d);
                                             
%% solve inverse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 1e-4; % regularization parameter

% construct forward solution map (matrix E)
source.Q(:) = 1; % this is to construct the forward map
[E, E2] = forwardErmak( source, recept, wind, ...
                        xmesh, ymesh, zmesh, dx, dy, dz, dt, tskip, nwind,...
                        Vspbo, Vdpbo, Mpb, stabclass, A, ...
                        IsAvgWind, IsOverDomain );
                  
% solve normal equation
Q = (E*E' + alpha*eye(source.n))\(E*d');


%% display solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('===============================================');
disp(['Solution         :', num2str(Q(1:6)')]);
disp([ num2str(Q(6:end)')]);
%disp(['Original         :', num2str(QQ(1:6)')]);
%disp([ num2str(QQ(6:end)')]);
disp('===============================================');
bar(Q*Mpb*365*24*60*60/1e3); % plot in tons/year

rmpath('.././readData/');
rmpath('.././data/');
rmpath('.././forwardsolve');
rmpath('.././data/inversedata');
rmpath('.././util');
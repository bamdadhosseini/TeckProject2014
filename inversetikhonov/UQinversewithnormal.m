%% Inverse solution using Tikhonov regularization
%% and Uncertainty quantification
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
                                             
%% solve inverse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 1e-4; % regularization parameter

% construct forward solution map (matrix E)
source.Q(:) = 1; % this is to construct the forward map
[E, E2] = forwardErmak( source, recept, wind, ...
                        xmesh, ymesh, zmesh, dx, dy, dz, dt, tskip, nwind,...
                        Vspbo, Vdpbo, Mpb, stabclass, A, ...
                        IsAvgWind, IsOverDomain );
            
nRealizations = 1000;                    
noiselevel= 0.4;
Q = zeros(source.n, nRealizations);                    
for ii = 1:nRealizations
   
    disp('realization:');
    disp(ii);
    % take a realization of data with noise
    % make sure the noise does not result in negative data
    while true

       variation = (rand(size(d))-1/2).*(noiselevel*(d));
       min(d + variation)
       if (d + variation > 0)
           dd(ii,:) = d + variation;
           break;
       end
    end
    
    % solve normal equation
    Q(:,ii) = (E*E' + alpha*eye(source.n))\(E*dd(ii,:)');
 
end

%% display solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qmean = mean(Q*Mpb*365*24*60*60/1e3,2);
Qvar = var(Q,1,2)*((Mpb*365*24*60*60/1e3)^2);
Qstd = std(Q*Mpb*365*24*60*60/1e3, 0, 2);

dmean = mean(dd,1);
dvar  = var(dd, 0, 1);
dstd = std(dd, 0  ,1);

disp('===============================================');
% disp(['Solution         :', num2str(Q(1:6)')]);
% disp([ num2str(Q(6:end)')]);
% %disp(['Original         :', num2str(QQ(1:6)')]);
% %disp([ num2str(QQ(6:end)')]);
disp('===============================================');
figure(1)
hb = bar(Qmean); % plot in tons/year
hold on
hv = errorbar(Qmean, sqrt(Qstd), 'xr');
set(gca, 'Xticklabel', source.label)
rotateXLabels(gca, 90);
ylabel('tn/yr');
%delete(hv(2)); % remove connecting lines

figure(2)
hd = bar(dmean); % plot in tons/year
hold on
hdv = errorbar(dmean, dstd, 'xr');
%delete(hv(2)); % remove connecting lines


rmpath('.././readData/');
rmpath('.././data/');
rmpath('.././forwardsolve');
rmpath('.././data/inversedata');
rmpath('.././util');
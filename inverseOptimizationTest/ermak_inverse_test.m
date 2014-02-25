%% ermak_inverse_test


clear all
clc
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
eps = 1e-3;
synthdep = totaldep + eps*rand(1,recept.n)*max(totaldep);
d = synthdep;
 
%%
%% Solve Inverse %%%%%%%%%%%%%%
%%%  uses steepest descent %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

alpha = 0;
Q = source.Qzn./max(source.Qzn);
Q = zeros(size(source.Qzn));
iter = 1;
error(1) = 2*optTol;
J(1) = (Q*E -d)*(Q*E -d)' + alpha*Q*Q';
QQ(1,:) = Q;
II = eye(size(E*E'));
errorQnorm = 2*optTol;
errorgradJ(1) = 2*optTol;
disp('Steepest Descent');
disp('');
gradJold = zeros(size(Q));

while (error(iter) >= optTol)
% while( error <= 20)   
    % Calculate Derivative
    %
    gradJ = (2*E*(Q*E- d)' + 2*alpha*Q')';
    
    
    % Choose search direction
    %
    Bk = eye(numel(gradJ));
    Pk = -1*gradJ*Bk';
    
    % Pick Optimal Step Size
    %
    beta = 1/2*(gradJ*gradJ')/(gradJ*(E*E' + alpha*II)*gradJ');
    % Update Solution
    %
    Q = Q + beta*Pk;
    iter = iter +1;
    
    % compute variation in solution of the 
    % functional
    J(iter) = (Q*E -d)*(Q*E -d)' + alpha*Q*Q';
    QQ(iter,:) = Q; 
    error(iter) = abs(J(iter) - J(iter-1))/abs(J(iter-1));
    %errorgradJ(iter) = norm( gradJold - gradJ, 2)/norm(gradJ,2);
    errorgradJ(iter) = norm(gradJ,2);
    gradJold = gradJ;
    errorQ = abs( QQ(iter,:) - QQ(iter-1,:) );
    disp(['iteration: ', num2str(iter), ' error: ', num2str(error(iter))]);
    disp(['iteration: ', num2str(iter), ' errorgradJ: ', num2str(errorgradJ(iter))]);
    %errorQnorm = norm(QQ(iter,:) - QQ(iter-1,:),2)/norm(QQ(iter-1,:),2)
end
% 
% QQ(1,:)
% QQ(end,:)
% synthsource.Qzn

disp('===============================================');
disp(['Original Source  :', num2str(synthsource.Qzn)]);
disp(['Solution         :', num2str(QQ(end,:))]);
disp(['Initial Guess    :', num2str(QQ(1,:))]);
disp('===============================================');

figure;
semilogy(2:iter, error(2:end), '-r');
xlabel('iteration');
ylabel('error');


figure;
semilogy(1:iter, J, '-b');
xlabel('iteration');
ylabel('J(Q)');

c= {'r', 'b', 'g', 'm'};
figure;
for i =1:source.n
    semilogy(1:iter, QQ(:,i), 'Color', c{i}  );
    hold on
end
xlabel('iteration');
ylabel('Q');
Legend('Q1', 'Q2', 'Q3', 'Q4');
hold off

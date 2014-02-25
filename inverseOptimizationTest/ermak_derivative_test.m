%% ermak_derivative_test

clear all
clc

%profile on
tic

setparams;
setSolver;

[depp, depp2] = forwardErmak(source, recept, wind, xmesh,...
    ymesh, zmesh, dx, dy, dz, dt, tskip, nwind, Vszns, Vdzns, Mzn, stabclass,A,...
    IsAvgWind, IsOverDomain);

totaldep = sum(depp,1); % sum over deposition from all sources to get
                       % the model output
data = totaldep + 0.0001*rand(1,recept.n);                      


%% compute derivative of functional %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J(Q) = \alpha*|| m - data ||^2_L2 + (1-\alpha)||Q||^2_L2 
% 
% where m is the output of model and alpha is the regularization
% parameter. Q is the solution which is the source emission rates.
% the derivative is then
% gradJ = 2*\alpha*(m - data)*[ d(vec(m))/d(vec(Q) ] + 2*(1 - \alpha)*Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute solution from each source with unit emission rate
% this will be used to construct the jacobian which in turn is 
% used in the derivative computation. This is infact the Jacobian
% [dm/dQ].
alpha = 1 - 0.1;

sourceUnit = source;
sourceUnit.Qzn = ones(size(source.Qzn));

[depUnit, depUnit2] = forwardErmak(sourceUnit, recept, wind, xmesh,...
    ymesh, zmesh, dx, dy, dz, dt, tskip, nwind, Vszns, Vdzns, Mzn, stabclass,A,...
    IsAvgWind, IsOverDomain);

% scale this result to compute total deposition in all receptors
dep = source.Qzn*depUnit;

% compute derivative

gradJ = 2*alpha*depUnit*(dep - data)' + 2*(1-alpha)*source.Qzn';

%% perform derivative test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we take random perturbations of the solution Q and evaluate
% the change in the functional. We reduce the size of perturbations 
% and expect convergence of at least order one for the derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the functional at current point
J = alpha*norm(totaldep - data, 2)^2 + (1-alpha)*norm(source.Qzn, 2)^2;

pertScale = 10.^[1:-1:-10]; % number of perturbations


for i = 1:numel(pertScale)
    pert = pertScale(i)*(rand(size(source.Qzn))-1/2);
    
    % solve the forward problem with perturbed values
    pertSource = source;
    pertSource.Qzn = pertSource.Qzn + pert;
    [pertDep, pertDep2] = forwardErmak(pertSource, recept, wind, xmesh,...
    ymesh, zmesh, dx, dy, dz, dt, tskip, nwind, Vszns, Vdzns, Mzn, stabclass,A,...
    IsAvgWind, IsOverDomain);
    totalpertDep = sum(pertDep,1);
    
    % compute variation in the functional
    pertJ = alpha*norm(totalpertDep - data,2)^2 + (1-alpha)*norm(pertSource.Qzn,2)^2;
    
    approxDiff(i) = abs(pertJ - J);
    approxDiffDer(i) = abs(pertJ - J - pert*gradJ);
end

approxDiff
approxDiffDer


%%%%%%%% PLOT %%%%%%%%%%%%%%
figure;
semilogy(1:numel(pertScale), approxDiff, 1:numel(pertScale), approxDiffDer);
xlabel('perturbation size');
ylabel('error');
legend('|f(x+s) - f(x)|', '|f(x+s) - f(x) - df*s|');
%%%%%%% show results %%%%%%%%%%%%%%
%data
%gradJ



toc
%profile viewer
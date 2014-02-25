%% computeGradient
%
% This function packages the adjoint solves to 
% return the gradient

function gradJ = computeGradient( sensor, source, wind,...
                    Vs, Vd, tskip, stabclass, R, alpha, C, targetQ)

       
      %% solve adjoint problem
      %
      % setup the adjoint source and solve the adjoint problem
      adjsource = getAdjointSourceForSensor( sensor, size(C));
      adjsol = solveAdjointConcentrationAtReceptors( sensor,...
                    source, wind, Vs, Vd, stabclass, ...
                    tskip, adjsource);
                
      %% compute gradient
      %
      gradJ = getGradientTikhonov( adjsol, source, R, alpha, targetQ);

end
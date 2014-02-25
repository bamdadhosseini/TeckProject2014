%% ArmijoLineSearch

function [Qk, gammak, Jk, C, sensor, lind, uind] = armijoLineSearch(sensor, source, wind, Vspbo, Vdpbo,...
    tskip, stabclass, R, alpha, targetQ, descentDir, J, gradJ, c1, tau, gammak, Armijomaxiter, i, l, u ) 
     %% choose step size (Armijo)
    
    disp('')
    disp('>>> Line Search')
    sourcek = source;
    LSiter = 1; %line search iteration
    
    while true
       [Qk, lind, uind] = constraintProjection(source.Q + gammak*descentDir, l, u); % positivity projection
       sourcek.Q = Qk;
       [Jk,C, sensor] = computeObjective(sensor, sourcek, wind, ...
                    Vspbo, Vdpbo, tskip, stabclass, R, alpha, targetQ);
       
       disp(['>>> J(Q_', num2str(LSiter), ')= ', num2str(Jk),...
            ' \gamma_k=', num2str(gammak), '  ', ...
             num2str(J(i) + c1*gammak*(descentDir(:)')*gradJ(:))]);
       %
       % first Wolfe Condition
       if (Jk <= J(i) + c1*gammak*(descentDir(:)')*gradJ(:)), break; end;
       %
       % Reduce step size
       gammak = gammak*tau;
       LSiter = LSiter + 1;
       
        if (LSiter > Armijomaxiter), disp('Max line search iteration reached!'), return; end;
    end
    
    if (LSiter == 1), gammak = 2*gammak; end;
end
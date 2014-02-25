%% getDescentDir with constraints

function [direction] = getDescentDirConstraint( grad, method, R, iter, l, u, ...
    sensor, source, wind, Vspbo, Vdpbo,...
        tskip, stabclass, alpha, targetQ, J, c1, tau, gammak, Armijomaxiter)
            
    global verbosity;

    if strcmpi(method,'STEEPESTDESCENT')
        direction = -grad;
        normfactor = 1./sqrt(sum(direction.^2, 2));
        %size(normfactor)
        direction = direction.*repmat(normfactor,1, size(grad,2));
        %pause;
        
    elseif strcmpi(method, 'GPRN')
        
        %% GPRN
        % 
        % use steepest descent to minimize using gradient projection
        stddir = -grad;
        normfactor = 1./sqrt(sum(stddir.^2, 2));
        stddir = stddir.*repmat(normfactor,1, size(grad,2));
        
        % call armijo to minimize
        disp('>>> Inner iteration <<<');
        [Qk, gammak, Jk, C, sensor, lind, uind] = armijoLineSearch(sensor, source, wind, Vspbo, Vdpbo,...
        tskip, stabclass, R, alpha, targetQ, stddir, J, grad, c1, tau, gammak, Armijomaxiter,iter, l, u );
        
        source.Q = Qk;
        %pause
        % get projected gradient at this location and find active sets
        projgrad = computeGradient(sensor, source, wind, ...
                    Vspbo, Vdpbo, tskip, stabclass, R, alpha, C, targetQ);
        
         %disp(size(lind))
         %disp(size(uind))
         %pause       
         
         % compute descent direction using reduced hessians
         for isrc =1:source.n
            
             % construct reduced hessian
             indx = [lind((lind(:,1)==isrc), 2); uind((uind(:,1)==isrc), 2)];
             e = zeros(size(R,1), 1);
             e(indx) = 1;
             Dact = spdiags( e, 0, size(R,1), size(R,2) );
             Dinact = spdiags( e == 0, 0, size(R,1), size(R,2) ); 
             
             projgrad(isrc, [lind((lind(:,1)==isrc), 2); uind((uind(:,1)==isrc), 2)]) = 0;
        
             
             %size(Dact)
             %size(Dinact)
             %size(R)
             %size(projgrad(isrc,:)')
             %size(Dinact*R*Dinact + Dact)
             %pause
             
             %figure(4)
             %spy(Dinact*R*Dinact + Dact)
             %pause
             direction(isrc, :) = - ( (Dinact*R*Dinact + Dact)\(projgrad(isrc,:)'));     
         end
        
        %direction = -((R)\(grad'))';
        normfactor = 1./sqrt(sum(direction.^2, 2));
        %size(direction)
        %size(normfactor)
        %pause
        direction = direction.*repmat(normfactor,1, size(grad,2));
        %pause
    else
        disp('Optimization method not available!');
        pause;
    end
end
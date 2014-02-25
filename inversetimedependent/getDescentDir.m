%% getDescentDir
%
% compute the descent direction for optimization algorithm

function [direction, s, y] = getDescentDir( grad, method, R, iter,...
                s, y, gradold, solold, sol, m, aa )
            
    global verbosity;

    if strcmpi(method,'STEEPESTDESCENT')
        direction = -grad;
        normfactor = 1./sqrt(sum(direction.^2, 2));
        %size(normfactor)
        direction = direction.*repmat(normfactor,1, size(grad,2));
        %pause;
        
    elseif strcmpi(method, 'REGULARIZATION')
        direction = -((R)\(grad'))';
        normfactor = 1./sqrt(sum(direction.^2, 2));
        direction = direction.*repmat(normfactor,1, size(grad,2));
        %direction = min(normfactor).*direction;
        
    elseif strcmpi(method, 'LBFGS') % see page 177-179 of Nocedal
       
        if (iter <= m) % we are initializing so use the Hessian of the
                      % regularization for the first guess
            if verbosity ==1
                disp(' ');
                disp(['>>>>initializing LBFGS : ', num2str(iter)]);
            end
            direction = -((R'*R)\(grad'))';
            normfactor = 1./sqrt(sum(direction.^2, 2));
            direction = direction.*repmat(normfactor,1, size(grad,2));
            
            s(:, iter) = sol(:) - solold(:);
            y(:, iter) = grad(:) - gradold(:);
%             disp(norm(s))
%             disp(norm(y))
%             pause;
        else           
            %norm(s)
            %norm(y)
            % use LBFGS
            q = grad(:);
            
            % first loop of LBFGS
            for ii = m:-1:1 
                rho(ii) = 1/(y(:,ii)'*s(:,ii));
                alpha(ii) = rho(ii)*(s(:,ii)'*q);
                q = q - alpha(ii)*y(:,ii);
            end
            
            % approximate hessian
            gammak = (s(:,end)'*y(:,end))/(y(:,end)'*y(:,end));
            %r = gammak*q;
            %disp(aa);
            r = ((R)\(grad'))';
            r = r(:);
            %r = (aa*(R'*R)\(grad'))';
            %r = r(:);
            %norm(r)
            %pause
            
            
            % second loop of LBFGS
            for jj = 1:m
               beta = rho(jj)*(y(:,jj)'*r);
               r = r + s(:,jj)*(alpha(jj)-beta);
            end
            % set direction
            direction = -reshape(r, size(sol,1), size(sol, 2));
            normfactor = 1./sqrt(sum(direction.^2, 2));
            direction = direction.*repmat(normfactor,1, size(grad,2));
            %direction = min(normfactor).*direction;
            % update memory
            s(:,1:end) = [s(:,2:end), sol(:) - solold(:)];
            y(:,1:end) = [y(:,2:end), grad(:) - gradold(:)];
%             norm(s)
%             norm(y)
%             alpha
%             rho
%             beta
%             
%             pause
        end
            
    else
        disp('Optimization method not available!');
        pause;
    end
end
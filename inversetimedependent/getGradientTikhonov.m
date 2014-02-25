%% getGradientTikhonov
%
% compute the gradient of the Tikhonov functional

function gradJ = getGradientTikhonov( adjsol, source, R, alpha, targetQ)

    global verbosity;
    
    QmT = source.Q - targetQ;
    gradR = QmT*R; % gradient of regularization term
    
    gradJ = adjsol + alpha*gradR;
    
    if (verbosity == 1) 
       
        disp(' ');
        disp(['||gradJ|| = ', num2str(norm(gradJ(:))),...
            ', ||gradR|| = ', num2str(norm(gradR(:)))]);  
        
    end
     
    
end
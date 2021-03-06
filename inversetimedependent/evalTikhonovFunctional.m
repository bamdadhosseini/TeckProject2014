%% evaluates the Tikhonov functional for current misfit and weights
%

function J = evalTikhonovFunctional( source, sensor, R, alpha, targetQ )
    
    global verbosity;

    aux = arrayfun( @(s) s.weight.*(sum(s.misfit.^2)), sensor, 'UniformOutput', 0);
    size(cell2mat(aux));
    aux = sum(cell2mat(aux));

    QmT = source.Q' - targetQ';
    
    RegTerm = R*QmT;
    RegTerm = reshape(QmT, numel(QmT), 1)'*RegTerm(:);
    %RegTerm = sum(sum(RegTerm.^2, 1));
%     max(max(RegTerm));
%     min(min(RegTerm));
    
    J = 0.5*aux + 0.5*alpha*RegTerm;
    
    if (verbosity ==1)
       
        disp(' ');
        disp(' ');
        disp(['||misfit||^2 = ', num2str(aux),...
            ', a*||R*Q||^2 = ', num2str(alpha*RegTerm)]);  
        
    end

end
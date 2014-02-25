% weighted Jacobi multi grid

function sol = weightedJacobiMgSolve( MG, init, RHS, nu1, nu2, w, tol, maxiter)

global verbosity;

sol = init;
iter = 1;
while (iter <= maxiter) || (res >= tol)
   
    [sol] = weightedJacobiMgIteration( MG, sol, RHS, nu1, nu2, w, 1);
    plot(sol)
    iter = iter +1;
    res = norm(MG(1).mat*sol - RHS)/length(RHS);  
    if verbosity == 1
       disp(['$$$ MG iteration : ', num2str(iter), '  ', num2str(res) ]);
    end
    pause
    
end


end
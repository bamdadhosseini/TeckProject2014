% single iteration of Jacobi Multi grid

function [sol] = weightedJacobiMgIteration( MG, uh, fh, nu1, nu2, w, level)

    sol = jacobiRelax(MG(level), uh, fh, nu1, w);
    plot(sol)
    pause
    if level == length(MG)

        sol = jacobiRelax(MG(level), uh, fh, nu2, w);
        return;
    else

        f2h = MG(level).If2c*( fh - MG(level).mat*uh);
        v2h = zeros(length(MG(level+1).mat),1);
        v2h = weightedJacobiMgIteration(MG, v2h, f2h, nu1, nu2, w, level+1);
        
    end
    
% correction
size(v2h)
size(MG(level).Ic2f)
max(uh)
pause
uh = uh + MG(level+1).Ic2f*v2h;
    
% relax again
sol = jacobiRelax(MG(level), uh, fh, nu2, w); 
    
    
end
% Jacobi Relaxation

%
% relax equation using Jacobi iterations

function sol = jacobiRelax( MG, u, f, nu, w )

    sol = u;
    size(sol)
    size(MG.RJ)
    for i =1:nu
        i
        %sol = ( (1-w)*speye(size(MG.RJ)) + w*MG.RJ )*sol + w*(MG.D\f);
       plot(sol, 'r')
        pause
        %sol = (MG.D)*( f + MG.RJ*sol );
        sol = MG.RG\(MG.U*sol + f);
        
    end

end
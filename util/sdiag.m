%% sdiag 
%
% construct a sparse diagonal matrix 

function smat = sdiag( vec )

    smat = spdiags( vec, 0, length(vec), length(vec));
    
end
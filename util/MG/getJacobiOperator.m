% get operators for Weighted Jacobi

function [D, RJ, RG, U] = getJacobiOperator( M )

% extract the operators from the matrix

D = diag(diag(M).^(-1));

LplusU = - triu(M,1) - tril(M,-1);

RJ = LplusU;

RG = (D + tril(M,1));
U  = -triu(M,-1);


end
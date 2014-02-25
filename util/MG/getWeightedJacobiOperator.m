% get operators for Weighted Jacobi

function [D, RJ] = getJacobiOperator( M )

% extract the operators from the matrix

D = spdiags(spdiags(M), 0, size(M,1), size(M,2));

LplusU = triu(M) + tril(M);

RJ = D\LplusU;

end
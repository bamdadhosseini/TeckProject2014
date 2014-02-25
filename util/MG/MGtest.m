close all
clear all
% Multi-Grid test script

N = 97;
dd = 1; % discretization size
R = defineRegularization(N, dd);

RR = R'*R; % this is the matrix we want to invert

% construct random solution

x = sin(linspace(0,2*pi,N)+ 3/2)';

% get RHS 
RHS = x;
%RHS = 2*ones(size(x));

%% now use MG to solve the system and retrieve x again

nMGlevels = 4; % number of MG levels 
nu1 = 6; % number of relaxation steps going up
nu2 = 6; %number of relaxations steps on residual coming up
tol = 1e-6;
maxiter = 500;


% get MG matrices (precompute for performance)
for i = 1:nMGlevels
    
    mat = defineRegularization(min(floor(N/(2^(i-1)))+1,N), i*dd); 
    MG(i).np  = floor((N-2)/(2^(i-1)))+2;
    %i
    %MG(i).np
    MG(i).mat = (mat*mat');
    MG(i).Ic2f = getMgOperator('coarse2fine', MG(i).np );
    MG(i).If2c = getMgOperator('fine2coarse', MG(i).np );
    [MG(i).D, MG(i).RJ, MG(i).RG, MG(i).U] = getJacobiOperator( MG(i).mat ); 
    
end

xxx = bicg(MG(1).mat, RHS, 1e-6, 10000);
plot(xxx, 'g');
pause
xx = jacobiRelax(MG(1), zeros(N,1), RHS, 10000, 0.2);
        plot(xx)
        figure
        plot(RHS)
        pause

w = 0.8;
xx = weightedJacobiMgSolve( MG, zeros(N,1), RHS, nu1, nu2, w, 1e-4, 10 );





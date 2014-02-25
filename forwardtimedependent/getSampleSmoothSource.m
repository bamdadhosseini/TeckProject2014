%% draw a sample from a smooth time varying random source

function Q = getSampleSmoothSource(nwind, sourcevariance)

% 
% Draws a sample from a smooth time varying source assuming 
% the prior knowledge that 
%
% Q^n = 1/2(Q^n-1 + Q^n+1) + normal(0, \sigma^2*I)
%
% input: 
%   nwind : number of time steps of the source (size of Q)
%   sourcevariance: \sigma for this source
%
% output:
%   Q : vector of source emission rate for each time step.
%

e = ones(nwind,1);
L = spdiags([-e, 2*e, -e], [-1,0,1], nwind, nwind);
Q = rand(1,1) +L\(sourcevariance*(randn(nwind,1)));
Q = Q.*(Q>=0);
end 
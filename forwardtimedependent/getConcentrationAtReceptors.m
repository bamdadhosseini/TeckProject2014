%% SolveForwardErmak
%
% solve the forward problem of finding concentrations at 
% each receptor in time using ermak solution and time dependent 
% source.
%
% input: 
%
%    - source : a structure containing information of sources as well as 
%               time dependent emission rates.
%    - recept : receptor information including their location.
%
%    - wind   : wind data, including time and everything.
%

function C = getConcentrationAtReceptors( sensor, source, ...
            wind, Vsettling, Vdeposition, stabclass, tskip)

% take time steps and compute the solution at each time 
nstep = [1:tskip:wind.n];
C = zeros(length(sensor), length(nstep));

% precompute distances and set up source and receptor heights (to improve
% performance)
%
% these matrices are used at each time step so we simply set them up 
% before they are used.
%

% compute distance of sources and receptors in each direction
%
% the matrix xdist_ij = x_recept(i) - x_source(j)
% similarly for ydist.

xdist = repmat([sensor(:).x]', 1, length(source.x)) - repmat(source.x', length(sensor), 1);
ydist = repmat([sensor(:).y]', 1, length(source.y)) - repmat(source.y', length(sensor), 1);

% rz and sz are matrices of repeated recept.z and source.z of size 
% [recept.n X source.n] that are used to compute the forward map.
%

rz    = repmat([sensor(:).z]', 1, source.n);
sz    = repmat(source.z', length(sensor), 1);

for k = nstep
   
    E = getForwardMapAtThisTime( xdist, ydist, rz, sz, wind.dir(k), wind.vel(k),...
        Vsettling, Vdeposition, stabclass );
    C(:, k) = E*source.Q(:,k);
    
end

end
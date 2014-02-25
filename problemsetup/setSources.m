%% function setSources
% sets up the source object with position and 
% label of sources ready for the solver

function [source] = setSources(pos, id)

source.n = length(pos.source);

source.x = (pos.source(:,1) - pos.origin(1,1)); % m
source.y = (pos.source(:,2) - pos.origin(1,2)); % m

% this is a compromise, we assume all sources are 
% at the same height since they represent an average
% of the sources in each area. This includes stacks and
% vents as well as roads and piles.
source.z = 5*ones(1,source.n)'; % m
source.label = id.source';     

% we have no estimate on the source emissions, consider 
% random vector.
tpy2kgps = 1.0 / 31536; 
source.Q = 20*rand(1,source.n)*tpy2kgps; % kg/s

end
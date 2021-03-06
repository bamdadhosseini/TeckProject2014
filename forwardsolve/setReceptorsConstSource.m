%% function setSources
% sets up the source object with position and 
% label of sources ready for the solver

function [recept] = setReceptorsConstSource(pos, id)

recept.n = length(pos.dustfall);

recept.x = (pos.dustfall(:,1) - pos.origin(1,1)); % m
recept.y = (pos.dustfall(:,2) - pos.origin(1,2)); % m

% this is a compromise, we assume all sources are 
% at the same height since they represent an average
% of the sources in each area. This includes stacks and
% vents as well as roads and piles.
recept.z = zeros(1,recept.n)'; % m (all are at ground)
recept.label = id.dustfall';     

end

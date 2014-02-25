%% function setReceptors
% sets the receptor object, this is not crucial but it is 
% needed as an emulator to use the old solver with the new 
% data structure.

function [recept] = setReceptors(sensor)

recept.n = length(sensor);

recept.x = [sensor(:).x]'; % m
recept.y = [sensor(:).y]'; % m

% this is a compromise, we assume all sources are 
% at the same height since they represent an average
% of the sources in each area. This includes stacks and
% vents as well as roads and piles.
recept.z = zeros(1,recept.n)'; % m (all are at ground)     

end
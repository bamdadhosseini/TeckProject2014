%% clearReceptors

% some of the receptors have problem during measurements. We leave these 
% out.

function newrecept = clearReceptors( recept, d )

cindx = find( d ~= 0);
newrecept.n = length(cindx);
newrecept.x = recept.x(cindx);
newrecept.y = recept.y(cindx);
newrecept.z = recept.z(cindx);
newrecept.label = recept.label(cindx);

end
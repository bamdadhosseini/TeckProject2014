%% positiveProjection

% apply a projection operator to satisfy positivity constraint

function v = positiveProjection( v )

v = v.*(v>=0);

end
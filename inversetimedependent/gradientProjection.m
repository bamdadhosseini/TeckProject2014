%% positiveProjection

% apply a projection operator to satisfy positivity constraint

function [v, lind, uind] = constraintProjection( v, l, u )

[lindr, lindc] = find(v <= l);
lind = [lindr, lindc];
v = v.*(v>=l) + l.*(v < l);

[uindr, uindc] = find(v >= u);
uind = [uindc, uindr];
v = v.*(v<=u) + u.*(v > u);

end
%% setXact

function Xact = setXact(pos, id)

Xact.n = length(pos.Xact);
Xact.x = (pos.Xact(:,1) - pos.origin(1,1)); % m
Xact.y = (pos.Xact(:,2) - pos.origin(1,2));
Xact.z = zeros(1, Xact.n)';
Xact.label = id.Xact';     

end
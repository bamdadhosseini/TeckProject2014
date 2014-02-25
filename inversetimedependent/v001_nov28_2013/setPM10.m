%% setPM10

function PM10 = setPM10(pos, id)

PM10.n = length(pos.PM10);
PM10.x = (pos.PM10(:,1) - pos.PM10(1,1)); % m
PM10.y = (pos.PM10(:,2) - pos.PM10(1,2));
PM10.z = zeros(1, PM10.n)';
PM10.label = id.PM10';     

end
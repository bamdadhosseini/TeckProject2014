%% setTSP

function TSP = setTSP(pos, id)

    for i=1:size(pos.TSP,1)
        TSP(i).x = (pos.TSP(i,1) - pos.origin(1,1)); % m
        TSP(i).y = (pos.TSP(i,2) - pos.origin(1,2));
        TSP(i).z = 0;
        TSP(i).label = id.TSP(i);     
    end

end
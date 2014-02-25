%% getMisfitAtSensor
%
% given the measurements and data at the sensor it computes the 
% data misfit.
%

function sensor = getMisfitAtSensor( sensor )

    % misfit is simply the difference between measurements and data 
    misfit = arrayfun( @(s) s.measurement - s.data, ...
        sensor, 'UniformOutput', 0);
    % distribute the computed misfit to corresponding sensor
    [sensor(:).misfit] = deal(misfit{:});

end
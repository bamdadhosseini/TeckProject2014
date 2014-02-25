%% getAverageConcentation( sensor, M, dt )


function measurement = getAverageConcentration( sensor )

    % loop over sensors
    for iSens = 1:length(sensor)

        % compute the average for each sensor based on its accumulation period
        measurement(iSens,:) = (1./[sensor(iSens).accum_period]').*...
            [sensor(iSens).accum_concentration];
    end
    
end
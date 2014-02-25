%% getSyntheticData( source, sensor, wind, Vs, Vd, stabclass, noiselevel )
%
% constructs synthetic data by generating random smooth sources
% and adding noise

function [sensor, source] = getSyntheticData( source, sensor, wind, Vs, Vd,...
                        stabclass, tskip, noiselevel )
   
    % construct a random matrix for source emissions (no real need for a smooth
    % one but it won't hurt)

    source.Q = zeros(source.n, wind.n);
    sourcevariance = 1e-3;
    for i =1:source.n
        source.Q(i,:) = getSampleSmoothSource(wind.n, sourcevariance); % mol/s 
    end
    %source.Q = source.Q + 1;
    figure(1);
    plot(wind.time(1:wind.n), source.Q);
    xlabel('time (matlab numeric)');
    ylabel('mol/s');     
    
    % solve forward problem to get data 
    C = getConcentrationAtReceptors(sensor, source, wind, Vs, Vd,...
                                    stabclass, tskip);
    sensor = getAccumulatedConcentration( C, sensor, 1:length(sensor));
    sensor = getMeasurementAtSensor( sensor );
    
    data   = arrayfun( @(s) abs(s.measurement + ...
        noiselevel*(mean(s.measurement))*(rand(size(s.measurement))-1/2)),...
        sensor, 'UniformOutput', 0);
    
    % distribute data to sensors 
    [sensor(:).data] = deal( data{:} );
end
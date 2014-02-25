%% GetMeasurementAtSensor
%
% take measurments of the accumulated concentrations at the sensors
%
% input  : sensor
%
% output : sensor
%

function sensor = getMeasurementAtSensor( sensor )

    % the sensor.scale is the constant which should be multiplied by 
    % accumulated concentrations in mol/s to get the proper measurements 
    % in the right units.
    meas = arrayfun( @(s) s.scale.*s.accum_concentration,...
                    sensor,'UniformOutput', 0);
                
    % I still have not figured out a way to skip the use of the 
    % auxiliary variable meas here. In order to distribute the 
    % computed values properly we need to use deal() which in turn
    % requires us to open up the cell array meas.
    
    [sensor(:).measurement] = deal(meas{:});  
    
end
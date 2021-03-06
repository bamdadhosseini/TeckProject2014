%% setSchedules
% 

% this function returns the schedule structure which determines the 
% schedule of each sensor including dutfall jars, Xacts, PM10s and TSPs 
% using the dates from the data file. 

function [schedule] = setSchedules( fnames, wind )

    % fnames = a cell array of the name of the variables to be loaded 
    % schedule = a structure including start and enc times as well as
    % accum_period of each sensor 
    
    for fn = fnames 
       
        data = load(fn);
        
        % figure out where the data belongs
        indx = find( sensor(:).label == data.label);
        
        % if there is more than one match there is a problem so cry for
        % help
        if numel(indx) >1; disp('Sensor label corresponds to more than one sensor');break;end;
        
        % for this data, figure out the starting times by comparison to
        % wind times 
        
        [time, start_indx, dummy_indx] = intersect( wind.time, data.time, 'row' );      
        
        schedule(indx).start_indx = start_indx;
    end

end
%% setSensorSchedule
%

% this functions sets the start and finish times of each sensor 
% measurement as well as their possible averaging periods 

function [sensor, index] = setSensorSchedule( sensor, wind, schedule, trustweight )

% identify sensors
% figure out index of different sensors based on their kind.
index.dustfall  = find(strcmp( {sensor(:).kind}, 'DUSTFALLJAR' ) );
index.xact     = find(strcmp( {sensor(:).kind}, 'XACT' ) );
index.TSP      = find(strcmp( {sensor(:).kind}, 'TSP' ) );
index.PM10      = find(strcmp( {sensor(:).kind}, 'PM10' ) );

% start with dustfall jars 
for i = index.dustfall

    % dustfall jars start measuring at the initial time and finish at the 
    % end of the period (measurements are for one month)
    
   sensor(i).start_indx   = 1;
   sensor(i).end_indx     = wind.n; 
   sensor(i).accum_period = wind.n; % dustfall jars do not take averages they only 
                                     % accumulate throughout the
                                     % simulation. We set these manually
                                     
end

% Xact 
for j = index.xact
   
   sensor(j).start_indx = 1;
   sensor(j).end_indx   = wind.n;
   sensor(j).accum_period = 1; % Xact takes hourly averages
end


% TSP
for k = index.TSP
   
    sensor(k).start_indx   = [1:10:700];
    sensor(k).end_indx     = [1:10:700];
    sensor(k).accum_period = 1; % TSP takes one hour averages as well
                                % (not sure about this but it will do 
                                % for these tests)
end

% PM10
for k = index.PM10
   
    sensor(k).start_indx   = [1:10:700];
    sensor(k).end_indx     = [1:10:700];
    sensor(k).accum_period = 1; % TSP takes one hour averages as well
                                % (not sure about this but it will do 
                                % for these tests)
end

trustweight(index.dustfall) =1;
trustweight(index.xact) = 1;
trustweight(index.TSP) = 1;
trustweight(index.PM10) = 1;
trustweight = num2cell(trustweight);


[sensor(:).weight] = deal(trustweight{:});

end
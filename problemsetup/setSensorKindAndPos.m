%% setSensors

%
% this is a more general code, here we consider all dustfalljars, Xact
% instruments, TSPs and PM10s as sensors. We assume these are devices which
% take accumulated concentration and somehow change it into data. In other
% words they all measure data of the type 
%           const*\int_t C(t)dt
%
% and the const determines whether it is accumulated deposition, 10minute 
% averages or anything else.
%



function sensor = setSensorKindAndPos( pos, id )

% input : 
%       -pos :vector of positions 
%       -is  :id of each sensor

% output :
%       -sensor : an array of structures each of which 
%                 contains information about the 
%                 position, scaling factor and start and finish times 
%                 of the sensors.
%

% start with the dustfall jars

for i=1:size(pos.dustfall,1)
 
    sensor(i).kind = ['DUSTFALLJAR'];
    sensor(i).x = pos.dustfall(i,1) - pos.origin(1,1);
    sensor(i).y = pos.dustfall(i,2) - pos.origin(1,2);
    sensor(i).z = 0;
    sensor(i).label = id.dustfall(i);
    
end

% append Xact sensors
n = length(sensor);
for i=1:size(pos.Xact,1)
 
    sensor(n + i).kind = ['XACT'];
    sensor(n + i).x = pos.Xact(i,1) - pos.origin(1,1);
    sensor(n + i).y = pos.Xact(i,2) - pos.origin(1,2);
    sensor(n + i).z = 0;
    sensor(n + i).label = id.Xact(i);
    
end

% append TSPs 
n = length(sensor);
for i=1:size(pos.TSP,1)
 
    sensor(n + i).kind = ['TSP'];
    sensor(n + i).x = pos.TSP(i,1) - pos.origin(1,1);
    sensor(n + i).y = pos.TSP(i,2) - pos.origin(1,2);
    sensor(n + i).z = 0;
    sensor(n + i).label = id.TSP(i);
    
end

%append PM10s
n = length(sensor);
for i = 1:size(pos.PM10,1)
 
    sensor(n + i).kind = ['PM10'];
    sensor(n + i).x = pos.PM10(i,1) - pos.origin(1,1);
    sensor(n + i).y = pos.PM10(i,2) - pos.origin(1,2);
    sensor(n + i).z = 0;
    sensor(n + i).label = id.PM10(i);    
end






end
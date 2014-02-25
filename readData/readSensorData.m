%% read excel data file and construct data structure

function data = readSensorData( filename )

    d = xlsread( filename);
    data.measurement = d(:,2);
    data.time = datevec(d(:,1));

end
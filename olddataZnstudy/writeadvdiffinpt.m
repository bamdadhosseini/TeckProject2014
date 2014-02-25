%% Script for reading wind data and writing as txt

% reads in data from .xls file and writes three seperate .data
% files containing the velocity, angle and time.

% specify .xls file address

%fname = '/home/visitor/bhosseini/AtmosProject/ConvCode/thrd0/data/WindDataJun3toJul2_2002B.xls'

%% read data into wind structure
% wind.vel : wind velocity in m/s
% wind.dir : wind direction measured from north in radians
% wind.time: time at which wind is measured

%wind = readwind(fname);
%wind.dir = wind.dir + ((170.+ 90)/180.*pi); %computational domain is rotated

%% write data as text

% write wind velocity
save('windvel.dat', '-struct', 'wind', 'vel', '-ascii', '-double');

% write wind direction
save('winddir.dat', '-struct', 'wind', 'dir', '-ascii', '-double');

% write time
save('windtime.dat', '-struct', 'wind', 'time', '-ascii', '-double');

%% END

%% read wind time, direction and speed

% All data was packaged in a single excel file. This is not 
% easy to deal with on all machines so we separate the data 
% on multiple .xls files which are easier to read into matlab.
% all files are stored under './data/'

function [wind] = readWind(ftimename, fveldirname)

%% wind time 
%
% read wind time first
time = xlsread(ftimename);

% convert to numeric values 
time = datenum(time);
% convert to date vector [year, month, day, hour, min, sec]
wind.time = datevec(time);

%% wind vel and dir
%
% read wind vel and dir
windveldir = xlsread(fveldirname);
wind.vel = windveldir(:,1);
wind.dir = windveldir(:,2);

%% save as a mat file
%
save( './data/windData', 'wind' ); 

end
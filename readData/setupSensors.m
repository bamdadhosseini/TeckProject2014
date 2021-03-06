%% Read data from PM10s, TSPs and Xact 
%
% this is a script, not a function. This is run once to read in the 
% data into matlab and then save it as a .mat file. This is useful 
% because all versions of MATLAB do not support newer excel files. we 
% only need to run this once with every new data


clear all
clc
close all
 
%addpath('../../readData');
addpath('../data/inversedata');

%% Birchbak %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PM10
BBpm10 =readSensorData('Birch_PM10_data.xlsx');
BBpm10.label = 'Birchbank PM10';

% TSP 
BBtsp =readSensorData('Birch_TSP_data.xlsx');
BBtsp.label = 'Birchbank TSP';

%% Butlerpart %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PM10
BPpm10 =readSensorData('BP_PM10_data.xlsx');
BPpm10.label = 'Butler Park PM10';

% TSP 
BPtsp =readSensorData('BP_TSP_data.xlsx');
BPtsp.label = 'Butler Park TSP';

%% Warzone %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PM10
WARpm10 =readSensorData('WAR_PM10_data.xlsx');
WARpm10.label = 'Warfield PM10';

%% Xact North (DF) %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% XACT
DFxact = readSensorData('DF_Xact_data.xlsx');
DFxact.label = 'Xact North';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE data as .mat files 
%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('../data/inversedata');
save('BBpm10.mat');
save('BBTSP.mat');
save('BPpm10.mat');
save('BPTSP.mat');
save('WARpm10.mat');
save('DFXact.mat');


rmpath('../data/inversedata');
%rmpath('../../readData');
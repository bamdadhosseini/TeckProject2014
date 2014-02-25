function [depdata] = readdustfall( fname )
% READRECEPTOR: Read the deposition data for each receptor
%       
%   Input parameters:
%
%         fname - Excel file name
%
%   Output parameter: a structure named 'depdata' with fields
%
%         depdata = deposition data

% Read the Excel file (it may have to be rewritten to an 
% older Excel format in order that Matlab doesn't choke on it). 


depdata = xlsread(fname);

end
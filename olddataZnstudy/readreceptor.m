function depdata = readreceptor( fname )
% READRECEPTOR: Read the deposition data for each receptor
%       
%   Input parameters:
%
%         fname - Excel file name
%
%   Output parameter: a structure named 'depdata' with fields
%
%         depdata.month(m)- month of measurements
%         depdata.day(m)  - day of the month measurement was taken
%         depdata.Zn(m,r) - zinc for month m, receptor r 
%         depdata.Sr(m,r) - strontium for month m, receptor r 
%         depdata.Su(m,r) - sulphur for month m, receptor r 
%         depdata.Ca(m,r) - calcium for month m, receptor r 

% Read the Excel file (it may have to be rewritten to an 
% older Excel format in order that Matlab doesn't choke on it). 
IsAvgExtraJars = 0;
warning( 'OFF', 'MATLAB:xlsread:Mode' );
[num,txt,raw] = xlsread( fname );
warning( 'ON',  'MATLAB:xlsread:Mode' );
len = size(num,1);

% NaN's correspond to trace amounts, so simply replace them with zeros.
[ii,jj] = find(isnan(num));
num(ii,jj) = 0.0;

dt = 0;
m  = 0;
for i = 1 : len,

  % The first time we hit each new date, parse the date into month
  % and day.
  if dt ~= num(i,1),
    dt = num(i,1);
    dtstr = num2str( dt );
    m  = m + 1; 
    depdata.Nm = m;
    depdata.Nr(m) = 0;
    depdata.dtstr(m,:) = dtstr;
    depdata.year(m)  = str2num(dtstr(1:4));
    depdata.month(m) = str2num(dtstr(5:6));
    depdata.day(m)   = str2num(dtstr(7:8));
    % When measurements overlap two months we need to avoid multiple
    % entries for a single month: when a measurement is taken at the
    % beginning of the following month, then use the previous month #.
    if depdata.day(m) < 5,
      depdata.month(m) = depdata.month(m) - 1;
    end
  end

  % Store the specific data we need for each receptor (but only the
  % first nine receptors which I have on the map).  Data values
  % given in mg must be converted to kg.  Only the "soluble"
  % fractions are used here. 
  r = num(i,2);   % dustfall jar number
  if r < 12,
    depdata.Nr(m) = depdata.Nr(m) + 1;
    depdata.Zn(m,r) = num(i,11) * 1e-6;
    depdata.Sr(m,r) = num(i,5)  * 1e-6; 
    depdata.Su(m,r) = num(i,13) * 1e-6;
    depdata.Ca(m,r) = 0.0;
  % For the duplicate dustfall jars (labelled 1/1001 and 5/1005) ... 
  elseif r > 1000, 
    rr = r-1000;
    if IsAvgExtraJars,
      % ... average the two measurements ...
      depdata.Zn(m,rr) = 0.5 * (depdata.Zn(m,rr) + num(i,11) * 1e-6); 
      depdata.Sr(m,rr) = 0.5 * (depdata.Sr(m,rr) + num(i,5)  * 1e-6); 
      depdata.Su(m,rr) = 0.5 * (depdata.Su(m,rr) + num(i,13) * 1e-6);
    end
    
  end
end

% END.

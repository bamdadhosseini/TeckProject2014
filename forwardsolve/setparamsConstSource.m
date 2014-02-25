% SETPARAMS: Set up various physical parameters for the atmospheric
%    dispersion problem.

grav  = 9.8;      % gravitational acceleration (m/s^2)
eta   = 1.8e-5;   % dynamic viscosity of air (kg/m.s)

% Densities all taken from Wikipedia:
rhozn = 7140;     % density of zinc (kg/m^3)
rhosr = 2640;     % density of strontium (kg/m^3)
rhosu = 2070;     % density of sulphur, alpha-form (kg/m^3) 
rhozns= 3540;     % density of zinc sulphate, ZnSO4 (kg/m^3)
rhosrs= 3960;     % density of strontium sulphate, SrSO4 (kg/m^3)
rhopb = 11340;    % density of lead (kg/m^3)
rhopbo= 9380;     % density of lead oxide (kg/m^3)

% Molar masses from Wikipedia, rounded to three decimal places:
Mzn   =  65.4e-3; % molar mass of zinc (kg/mol)
Msr   =  87.6e-3; % molar mass of strontium (kg/mol)
Msu   =  32.1e-3; % molar mass of sulphur (kg/mol) 
Mso4  =  96.1e-3; % molar mass of sulphate ion (kg/mol)
Mzns  = 162.0e-3; % molar mass of zinc sulphate (kg/mol)
Msrs  = 184.0e-3; % molar mass of strontium sulphate (kg/mol)
% from convertunits.com/molarmass ( they have a nice calculator )
Mpb   = 207.2e-3; % molar mass of lead (kg/mol)
Mpbo  = 223.1e-3; % molar mass of lead oxide PbO(kg/mol)

dzn   = 0.9e-6;   % diameter of zinc particles (m) - see Gatz (1975) 
dsr   = dzn;      % diameter of strontium particles (m)
dsu   = dzn;      % diameter of sulphur particles (m)
dzns  = 5.0e-6;   % diameter of zinc sulphate particles (m) - see CuSO4 in Sehmel (1980)
dsrs  = 5.0e-6;   % diameter of strontium sulphate particles (m) - same as ZnSO4
% from teck report (see Environmental dustfall partifulate characterization
% report for Oct 18,2013)
dpbo = 3.0e-6; % mean of particle sizes for PbO in the measurements.

Vdzn = 0.0062;    % Zn deposition velocity (m/s) -- in the range [5e-4,1e-2] m/s.  
Vdsr = 0.030;     % See Gatz, McMahon & Denison, Chrysikopoulos &c, 
Vdsu = 0.008;     % Lin & Hildemann, Rasmussen &c.
                  % [ Sulphur is approximated using SO2, Sr by Ca.
		  %   Sehmel (1980) lists Zn in [0.004,0.045] m/s and 0.018
		  %   and Sr in [2e-5,1e-4]. ]
Vdzns= 0.005;     % deposition velocity for ZnSO4 (m/s) - see Pacyna &c (1989) and Sehmel.
Vdsrs= 0.005;     % deposition velocity for SrSO4 (m/s) - same as ZnSO4 
Vdpbo= 0.005;

% Settling velocities from Stokes' law (m/s):
Vszn = rhozn * grav * dzn^2 / 18 / eta; 
Vssr = rhosr * grav * dsr^2 / 18 / eta; 
Vssu = rhosu * grav * dsu^2 / 18 / eta;  
Vszns= rhozns* grav * dzns^2/ 18 / eta;  
Vssrs= rhosrs* grav * dsrs^2/ 18 / eta;  
Vspbo = rhopbo* grav * dpbo^2/18 / eta; 
				 
dr = 0.162;       % receptor container diameter (m)
A  = pi*(dr/2)^2; % receptor container area (m^2) 

% read position data and id of all sources and receptors.
addpath('./data');
TeckPosData;
rmpath('./data');

% set position and label of sources (Note these are 
% area sources approximated by a single point source 
% in the centre.
source = setSources(pos, id);

% set position and label of dustfall jars. There are 
% 26 dustfall jars set on the area. Note these measure
% deposition and are different from TSP and Xact real time
% monitors.
recept = setReceptorsConstSource(pos, id);

% set Xact TSP and PM10
%

%Xact = setXact( pos, id);
%PM10 = setPM10( pos, id);
%TSP  = setTSP ( pos, id);

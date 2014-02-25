%% Set problem 

isanim = 0;         % animate the plumes
isplot = 0;         % plot deposition (takes a LONG time if tskip=1!)
tskip =  1;        % skip time steps in deposition calculation
IsAvgWind = 0;      % average the wind data if tskip > 1?
IsBWplots = 1;      % write black/white plot versions?
mycolormap = summer;% suggestions: summer, jet, hot, gray, ...
mybwcolmap = gray;   


umin = 0.1;
nrec = 26;
stabclass = 'c';

% Set up other solver parameters. 
nx = 50;
ny = nx;
nz = 50;
xlim = [-2000,4000];
ylim = [-2000,4000];
zlim = [0,100];
dx = (xlim(2)-xlim(1)) / nx;
dy = (ylim(2)-ylim(1)) / ny;
dz = (zlim(2)-zlim(1)) / nz;
x0 = [xlim(1):dx:xlim(2)]; % distance along the wind direction (m)
y0 = [ylim(1):dy:ylim(2)]; % cross-wind distance (m)
z0 = [zlim(1):dz:zlim(2)]; % vertical distance (m)
[xmesh,ymesh] = meshgrid(x0,y0);
zmesh = 0;    
smallfont = [14];

pskip = 100;  % skip time steps for animations only


IsAvgWind = 0;   % average the wind data if tskip > 1?
IsOverDomain = 0; % do we want to have the solution over whole domain?
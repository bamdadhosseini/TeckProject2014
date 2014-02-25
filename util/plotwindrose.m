function plotwindrose( wind, umin )
% PLOTWINDROSE: Using the given wind information, plot the wind rose
%   including both direction and velocity data. 
%       

figure;

global vnbin vbin vedges;

nbar = 16;
vedges  = [0, umin+1e-5, 1.0, 2.0, 3.0, 4.0, inf];
nvedges = length(vedges);
umax = max(wind.vel);

% First, calculate the percentage of calms.
[vnbin, vbin] = histc( wind.vel, vedges );
ntot = length(wind.vel);
rcalm = vnbin(1)/ntot*100;

% Next, divide the angular data into bins and select the largest 
% one to determine the maximum radial extent.
dir2 = wind.dir  - pi/nbar;
thedges = ([0:nbar]-0.5)/nbar * (2*pi);
[thnbin, thbin] = histc( dir2, thedges );
rmax  = max(thnbin)/ntot*100;
r0 = 0.111*rmax;
rmax  = rmax + r0;
if rmax/5 >= 3,
  rmult = 10;
else
  rmult = 5;
end
rmax = floor(rmax/rmult)*rmult;

% Plot dotted circles at major divisions in percentages.
%figure(1)
clf
for i = 1 : ceil(rmax/rmult),
  plotcircle( 0, 0, i*rmult, 'k:' )
end

% Plot the NSEW arrows.  
hold on
plot( [-1.1*rmax,-r0,0,r0,1.1*rmax], [0,0,NaN,0,0], 'k-', 'Linewidth', 0.2 )
plot( [0,0,0,0,0], [-1.1*rmax,-r0,NaN,r0,1.1*rmax], 'k-', 'Linewidth', 0.2 )
hold off
text( -0.04*rmax,  1.15*rmax, 'N' )
text( -0.04*rmax, -1.20*rmax, 'S' )
text( -1.27*rmax, -0.03*rmax, 'W' )
text(  1.15*rmax, -0.03*rmax, 'E' )
text( 0.90*rmax, -0.1*rmax, [num2str(rmax),'%'], 'Fontsize', [16] );
text( 0.40*rmax,  0.9*rmax, ['Max. wind=', num2str(umax,4),' m/s'], 'Fontsize', [16] );
axis off
axis equal

% Plot a circle representing the "calm" wind measurements.
hold on
plotcircle( 0, 0, r0, 'k-' )
text( -0.08*rmax, -0.03*rmax, [num2str(rcalm,2),'%'], 'Fontsize', [14] )
hold off

for i = 1 : nbar,
  
  jj = find(thbin==i);
  vv = wind.vel(jj);
  [vnbin, vbin] = histc( vv, vedges );
  vnbin = vnbin / ntot * 100;  % a percentage of total wind values

  theta0 = (i-1)/nbar * (2*pi);

  pclast = r0;
  for j = 2 : (nvedges-1);
     r = pclast + [0, vnbin(j), vnbin(j), 0, 0];
     thick  = 0.75;      % thickness of bars in [0,1]
     afac   = 2.0/thick; % multiplier in [2,inf]
     theta1 = theta0 - (j-1)*(2*pi)/nbar/(nvedges-2)/afac;
     theta2 = theta0 + (j-1)*(2*pi)/nbar/(nvedges-2)/afac;
     theta = [theta1, theta1, theta2, theta2, theta1];
     x = r .* cos(theta);
     y = r .* sin(theta);
     hold on
     ihsv = [ 0.667, j/(nvedges-1), j/(nvedges-1)];
     irgb = rgb2hsv(ihsv);
     colormap jet;
     aa = [ 1.0, 0.0, 0.0;
	    0.8, 0.0, 0.2;
	    0.6, 0.0, 0.4;
	    0.4, 0.0, 0.6;
	    0.2, 0.0, 0.8;
	    0.0, 0.0, 1.0
	  ];
     aa=colormap;
     irgb = aa( 65-floor(j/(nvedges-1)*36), :);
     %%irgb = aa( j, :);
     hvec(j-1) = fill( x, y, irgb);
     hold off
     
     pclast = pclast + vnbin(j);
  end

end

legend boxoff
% uminstr = sprintf( '%1d', floor(umin) ); % cut off decimal for umin
uminstr = sprintf( '%3.1f', umin );
legend( hvec, [ [' ',uminstr,'-1']; ...
		'  1-3 '; '  2-3 '; ...
		'  3-4 '; '  >4  ' ], ...
	'Location', 'SouthWest' );
hold off

% END PLOTWINDROSE.



function plotcircle( xc, yc, r, ltype )

if nargin < 4,
  ltype = 'k:';  % default line type
end

npts = 200;
theta = [0 : 2*pi/npts : 2*pi];
hold on
plot( xc+r*cos(theta), yc+r*sin(theta), ltype, 'Linewidth', 0.2 )
hold off

% END PLOTCIRCLE.

function C = ermakvec( x, y, z, H, Vw, Vs, Vd, stabclass )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is a vectorized version of John's ermak solver 
% to construct forward maps at each time step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ERMAK: Compute the concentration of pollutant (in mol/m^3) using the
%   Gaussian plume model modified for a deposition and settling
%   velocity.  This code handles one source (Q, located at the
%   origin) and multiple receptors (with locations given by x,y,z). 
%
% Input parameters:
%
%     x - receptor location: distance along the wind direction, with
%         the source at x=0 (m) 
%     y - receptor location: cross-wind direction (m)
%     z - receptor location: vertical distance (m)
%     H - height of source (m)
%     Q - pollutant emission rate (mol/s)
%    Vw - wind velocity (m/s)
%    Vs - gravitational settling velocity (m/s)
%    Vd - deposition velocity (m/s)
% stabclass - stability class (A-F, optional)   
%
% Output:
%
%     C - contaminant concentration (mol/m^3)
%
% References: Ermak (1977), Winges (1990/1992).

% First, define the velocity cut-off, below which the concentration
% is set to zero.
Vcutoff = 0.0;
sigmatype = 2;

% Next are Brookhaven numbers for stability classes.  
% Source: Canepa (coefficients), epa.gov (class vs wind, based on daytime).
%
% Pasquill-Gifford coefficients:  sigma = a*|x|^b
% Briggs coefficients:            sigma = c*x / (1+d*x)^e
%
% Source: Carrascal, Puigcerver & Puig
%         (Theor. Appl. Climatol. 48:147-157, 1993).  

% Coefficients are based on stability class. 
%
% The default class is the "neutral" stability class D (Ed's
% suggestion).  This is a better approximation for low wind conditions
% where there is no inversion and they see "particulate hanging".  This
% is a good approximation except for the month(s) ending in November.
if nargin < 9,
  stabclass = 'D';
end
switch lower(stabclass)
  case 'a'
   ay = 0.40;  by = 0.91;   az = 0.41; bz = 0.91;  % A: very unstable (0-2 m/s)
   cy = 0.22;  dy = 0.0001; ey = 0.5;  
   cz = 0.20;  dz = 0.0;    ez = 0.0;
 case 'b'
  ay = 0.36;  by = 0.86;   az = 0.33; bz = 0.86;   % B: unstable (2-3 m/s)
  cy = 0.16;  dy = 0.0001; ey = 0.5; 
  cz = 0.12;  dz = 0.0;    ez = 0.0;
 case 'c'
  ay = 0.36;  by = 0.86;   az = 0.33; bz = 0.86;   % C: slightly unstable (3-5 m/s)
  cy = 0.11;  dy = 0.0001; ey = 0.5; 
  cz = 0.08;  dz = 0.0002; ez = 0.5;
 case 'd'
  ay = 0.32;  by = 0.78;   az = 0.22; bz = 0.78;   % D: neutral (5+ m/s)
  cy = 0.08;  dy = 0.0001; ey = 0.5; 
  cz = 0.06;  dz = 0.0015; ez = 0.5;
 case 'e'
  ay = 0.315; by = 0.745;  az = 0.14; bz = 0.745;  % E: slightly stable
  cy = 0.06;  dy = 0.0001; ey = 0.5;
  cz = 0.03;  dz = 0.0003; ez = 1.0;
 case 'f'
  ay = 0.31;  by = 0.71;   az = 0.06; bz = 0.71;   % F: stable
  cy = 0.04;  dy = 0.0001; ey = 0.5;
  cz = 0.016; dz = 0.0003; ez = 1.0;
 otherwise
  %%ay = 0.31;  by = 0.71;   az = 0.06; bz = 0.71;   % G: ?
  error('Invalid stability class (should be A-F)')
end

% ay = 0.32;  by = 0.78;   az = 0.22; bz = 0.78;   % D: neutral (5+ m/s)
% cy = 0.08;  dy = 0.0001; ey = 0.5; 
% cz = 0.06;  dz = 0.0015; ez = 0.5;

if sigmatype == 1,
  % Pasquill-Gifford:
  sigmay = ay*abs(x).^by .* (x > 0);
  sigmaz = az*abs(x).^bz .* (x > 0);
else
  % Briggs:
  sigmay = cy*abs(x) .* (1.0 + dy*x).^(-ey) .* (x > 0);
  sigmaz = cz*abs(x) .* (1.0 + dz*x).^(-ez) .* (x > 0);
end

% Calculate the eddy diffusivity based on the assumptions 
% K is the same in all directions and is approximately constant
% (Winges, p. 5):  
K = 0.5*sigmaz.^2.*Vw./x;  % m^2/s

% Calculate the concentration of pollutant (in mol/m^3).
if Vw < Vcutoff,
  C = 0 * x;
else
  V1 = Vd - 0.5*Vs;
  C  = 1 ./ (2*pi*Vw*sigmay.*sigmaz) .* exp( -0.5*y.^2./sigmay.^2 ) .* ...
       exp( -0.5*Vs*(z-H)./K - Vs^2*sigmaz.^2/8./K.^2 ) .* ...
       ( exp( -0.5*(z-H).^2./sigmaz.^2 ) + ...
	 exp( -0.5*(z+H).^2./sigmaz.^2 ) - sqrt(2*pi)*V1*sigmaz./K .* ...
	 exp( V1*(z+H)./K + 0.5*V1^2*sigmaz.^2./K.^2 ) .* ...
	 erfc( V1*sigmaz/sqrt(2)./K + (z+H)./sqrt(2)./sigmaz ) );
  ii = find(isnan(C) | isinf(C));
  C(ii) = 0;   % Set all NaN values to zero.
end

% END.
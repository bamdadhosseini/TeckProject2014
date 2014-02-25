%% getDepDomain
%
% compute deposition at receptors

function dep2 = getDepDomain( dep2, xmesh, ymesh, zmesh, dx, dy,...
     source, theta, Vw, Vs, Vd, stabclass,...
     M, dt)

 % getDepAtREcept - compute deposition at receptors only
 %
 % input :-
 %
 % xmesh     - x comp of meshgrid
 % ymesh     - y comp of meshgrid
 % dx        - increment in x
 % dy        - increment in y
 % source    - structure containing source location and emissions
 % theta     - angle of rotation of the ermak solutions
 % Vw        - wind velocity in ermak solution
 % Vs        - settling velocity
 % Vd        - deposition velocity
 % stabclass - atmospheric stability class
 % M         - Molar mass of pollutant
 % dt        - timesteps in model
 
 dep2i = 0;
 for i=1:source.n

      xxx = (xmesh - source.x(i))*cos(theta) + (ymesh - source.y(i))*sin(theta);
      yyy =-(xmesh - source.x(i))*sin(theta) + (ymesh - source.y(i))*cos(theta);
      zzz = zmesh;
      C2 = ermak( xxx, yyy, zzz, source.z(i), source.Q(i), ...
                   Vw, Vs, Vd, stabclass );
      % C2 is a concentration in mol/m^3, convert to deposition in kg ...
      dep2i = dep2i + ((dx*dy) * dt * M * Vd) * C2; 
    
 end
 dep2(:,:, i) = dep2(:,:,i) + dep2i;
end
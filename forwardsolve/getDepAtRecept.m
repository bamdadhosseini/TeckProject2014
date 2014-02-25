%% getDepAtRecept
%
% compute deposition at receptors

function dep = getDepAtRecept( dep, recept, source, theta, Vw, Vs, Vd, stabclass,...
     Arec, M, dt)

 % getDepAtREcept - compute deposition at receptors only
 %
 % input :-
 %
 % recept    - structure containing receptor locations
 % source    - structure containing source location and emissions
 % theta     - angle of rotation of the ermak solutions
 % Vw        - wind velocity in ermak solution
 % Vs        - settling velocity
 % Vd        - deposition velocity
 % stabclass - atmospheric stability class
 % Arec      - area of receptor (dustfall jar) 
 % M         - Molar mass of pollutant
 % dt        - timesteps in model
 
 
 for i=1:source.n
    % First, accumulate depositions at the receptor sites.
    xxx = (recept.x - source.x(i))*cos(theta) + (recept.y - source.y(i))*sin(theta);
    yyy =-(recept.x - source.x(i))*sin(theta) + (recept.y - source.y(i))*cos(theta);
    C = ermak( xxx, yyy, recept.z, source.z(i), source.Q(i), ...
               Vw, Vs, Vd, stabclass );
    % C is a concentration in mol/m^3, convert to deposition in kg ...
    %size(C)
    %size(dep)
    %pause
    dep(i,:) = dep(i,:) + (Arec * dt * M * Vd) * C'; 

 end
end

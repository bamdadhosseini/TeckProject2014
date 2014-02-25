%% function for solution of the forward problem

function [dep, dep2] = forwardErmak(source, recept, wind, ...
                            xmesh, ymesh, zmesh, dx, dy, dz, dt, tskip, nwind,...
                            Vs, Vd, M, stabclass,A,...
                            IsAvgWind, IsOverDomain )

    dep  = zeros(source.n, recept.n);
    dep2 = zeros(size(xmesh,1), size(xmesh,2), source.n);% zeros(source.n,size(xmesh));
    
%     for k = [1 : tskip : nwind],
% 
%       k2 = min( k+tskip-1, nwind );
%       if IsAvgWind,
%         thetaavg = mean( wind.dir(k:k2) );
%         Vwavg    = mean( wind.vel(k:k2) );
%       else
%         thetaavg = wind.dir(k);
%         Vwavg    = wind.vel(k);
%       end
%       
%       %if floor((k-tskip)/500) ~= floor(k/500), fprintf( 1, ' %5d ..', k ), end
%       theta = +(pi/2 - thetaavg);    % ** CHECK THIS! **
%       Vw    = Vwavg;
%       warning( 'OFF', 'MATLAB:divideByZero' );
% 
%       for i = 1 : source.n,
%         %% First, accumulate depositions at the receptor sites.
%         xxx = (recept.x - source.x(i))*cos(theta) + (recept.y - source.y(i))*sin(theta);
%         yyy =-(recept.x - source.x(i))*sin(theta) + (recept.y - source.y(i))*cos(theta);
%         C = ermak( xxx, yyy, recept.z, source.z(i), source.Qzn(i), ...
%                    Vw, Vszns, Vdzns, stabclass );
%         % C is a concentration in mol/m^3, convert to deposition in kg ...
%         dep(i,:) = dep(i,:) + (A * dt * Mzn * Vdzns) * C; 
% 
%         if IsOverDomain,
%           %% Next, accumulate depositions in all cells on the mesh.
%           xxx = (xmesh - source.x(i))*cos(theta) + (ymesh - source.y(i))*sin(theta);
%           yyy =-(xmesh - source.x(i))*sin(theta) + (ymesh - source.y(i))*cos(theta);
%           zzz = zmesh;
%           C2 = ermak( xxx, yyy, zzz, source.z(i), source.Qzn(i), ...
%                        Vw, Vszns, Vdzns, stabclass );
%           % C2 is a concentration in mol/m^3, convert to deposition in kg
%           % ...
%           dep2i = ((dx*dy) * dt * Mzn * Vdzns) * C2; 
%         end
%       dep2(:,:, i) = dep2(:,:, i) + dep2i;
%       %dep2i = 0;
%       end
%       warning( 'ON', 'MATLAB:divideByZero' );
%     end
%     
    for k = [1 : tskip : nwind],

      k2 = min( k+tskip-1, nwind );
      if IsAvgWind,
        thetaavg = mean( wind.dir(k:k2) );
        Vwavg    = mean( wind.vel(k:k2) );
      else
        thetaavg = wind.dir(k);
        Vwavg    = wind.vel(k);
      end

      %if floor((k-tskip)/500) ~= floor(k/500), fprintf( 1, ' %5d ..', k ), end
      theta = +(pi/2 - thetaavg);    % ** CHECK THIS! **
      Vw    = Vwavg;
      warning( 'OFF', 'MATLAB:divideByZero' );
      
      % calculate deposition at receptors
      dep = getDepAtRecept( dep,recept, source, theta, Vw, Vs, Vd, stabclass ...
                             ,A, M, dt);
      if IsOverDomain,
      % Next, accumulate depositions in all cells on the mesh.
        dep2 = getDepDomain( dep2, xmesh, ymesh, zmesh, dx, dy, source, theta,...
                                    Vw, Vs,Vd, stabclass, M, dt); 
      end
      warning( 'ON', 'MATLAB:divideByZero' );

    end


end
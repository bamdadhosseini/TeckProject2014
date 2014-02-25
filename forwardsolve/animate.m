%% script for animation.
%
% this is under developed. simply cut out of the main script to improve
% readablility.

  if isanim && (mod(k-1+pskip,pskip) == 0 || k == nwind), 
    theta = +(pi/2 - thetaavg);
    Vw    = Vwavg;

    warning( 'OFF', 'MATLAB:divideByZero' );
    C = 0;
    for i = 1 : source.n,
      xxx = (xmesh - source.x(i))*cos(theta) + (ymesh - source.y(i))*sin(theta);
      yyy =-(xmesh - source.x(i))*sin(theta) + (ymesh - source.y(i))*cos(theta);
      zzz = zmesh;
      C = C + ermak( xxx, yyy, zzz, source.z(i), source.Q(i), ...
                     Vw, Vspbo, Vdpbo, stabclass );   % in mol/m^3
    end

    warning( 'ON', 'MATLAB:divideByZero' );
    
    if max(C(:)) > 0,
      contour( xmesh, ymesh, C, ncontour );
      %contour3( xmesh, ymesh, C, 20 );
      set(gca,'xlim',[min(xmesh(:)),max(xmesh(:))]);
      set(gca,'ylim',[min(ymesh(:)),max(ymesh(:))]);
      xlabel('x (m)'), ylabel( 'y (m)' )
      colorbar
      title( sprintf( 'Concentration, mol/m^3 (t=%6.2f h)', wind.time(k)/3600 ) )
      grid on  shg
      pause(0.01)   % force each frame to appear.
    end
  end
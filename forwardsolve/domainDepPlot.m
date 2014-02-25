  
% plot deposition contours 


  figure(3), clf
   clist = [ 0.01, 0.04, 0.06, 0.08, 0.1, 0.015]; 
    [cdep2, hdep2] = contourf( xmesh, ymesh, dep2/(dx*dy)*1e3, clist );
    max(max(max(dep2/(dx*dy)*1e3)))
%   [cdep2, hdep2] = contourf( xmesh, ymesh, dep2/(dx*dy)*1e3, 10);
%    clabel(cdep2, hdep2, 'FontSize', smallfont-2 )
  %colormap jet;   % colormap(0.8*(1-jet))
  colormap(1-winter);  % for filled contours
  colorbar
  hold on
  plot( source.x, source.y, 'ro', 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'r' )
  text( source.x, source.y, source.label, 'FontSize', smallfont, 'FontWeight','bold' );
  axis equal
  set(gca, 'XLim', xlim ), set(gca, 'YLim', ylim )
  xlabel('x (m)'), ylabel('y (m)')
  title('Contours of Zn concentration (g/m^2)')
  greencolor = [0,0.816,0];  % matches triangles in trailsite3.eps
  plot( recept.x(1:nrec), recept.y(1:nrec), 'g^', 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', greencolor )
  text( recept.x(1:nrec), recept.y(1:nrec), recept.label(1:nrec), ...
        'FontSize', smallfont, 'FontWeight', 'bold' )
  hold off
  if IsBWplots,
    colormap(gray)
    print -depsc 'dep1bw.eps'
  end
  colormap(1-winter);  % for filled contours
  print -depsc 'dep1.eps'
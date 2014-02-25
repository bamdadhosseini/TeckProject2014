% INVERSE: Solve the inverse atmospheric dispersion problem.  
%    That is, given the amount of Zn deposited in each receptor
%    (in kg), estimate the source emission rates (in mol/s). 

tic
clear all
nContam    = 3;     % 1 = Zn only, 2 = Zn+Sr, 3 = Zn+Sr+SO4
IsPause    = 0;
IsAllMonths= 0;     % for multiple months, optimize over all months simultaneously
umin       = 0.1    % wind cut-off (0.5 according to Hanna et al.(1982))
dnoise     = 0.0;   % add noise?
IsUseSr    = 1;     % include Sr in the equality constraints?
IsPlotMeanD= 0;     % add avg deposition (avgd over months 1,4:8) to bar plot
IsPlotZnEst= 0;     % plot a point for the engineering estimates of Zn sources
IsAvgWind  = 0;     % average the wind data if tskip > 1?
IsBWplots  = 1;     % write black/white plot versions?
stabclass  = 'D';   % stability classes A-F (D=neutral, E=slightly stable)
tskip      = 1;     % increase to 20 if a fast run is desired

mycolormap = summer;% suggestions: summer, jet, hot, gray, ...
mybwcolmap = gray;   

nr    = 9;          % number of receptors to include (9 or 11)
rlist = [1,2,4:nr]; %   - omit R3  **OR**
rlist = [1:nr];     %   - include all receptors

mlist = [1,5,6,7];  % Index of months in deposition data file --
                    % "**" refers to months where wind data is available:
                    % 1 = Jun 1-29, 2001      **
                    % 2 = Jun 30-Jul 26, 2001 [anomalously low values]
                    % 3 = Jul 27-Aug 30, 2001 [ "" ]          
                    % 4 = ???-Oct 24, 2001   
                    % 5 = Oct 25-Nov 23, 2001 **
                    % 6 = May 1-30, 2002      **
                    % 7 = Jun 3-Jul 2, 2002   **  (DEFAULT)
                    % 8 = Jul 3-30, 2002
                    % negative = take deposition equal to average over
                    %     all months (except Jul/Aug 2001) but wind 
                    %     data from this month.
mlist = [1,4:8];    % skip anomalously low measurements
mlist = [1,4,5];    % [1,4,5] for 2001, [6,7,8] for 2002
mlist = [6,7,8];    % [1,4,5] for 2001, [6,7,8] for 2002
mlist = 7;          % default: June 2002
mlist = [1:8];      % all months

setparams;
depdata = readreceptor( 'LIMS_Data_EW48.xls' );

%source.z = source.z * 1.4;
%recept.z = recept.z * 1.4;
%fac = 1.4; Vdzn=Vdzn*fac; Vdsu=Vdsu*fac; Vdsr=Vdsr*fac; Vdzns=Vdzns*fac; Vdsrs=Vdsrs*fac;

% Exact dates
monthstr = ['Jun 1-29, 2001     '; 'Jun 30-Jul 26,2001 '; ...
            'Jul 27-Aug 30, 2001'; 'Sep 27-Oct 24,2001 '; ...
            'Oct 25-Nov 23, 2001'; 'May 1-30, 2002     '; ...
            'Jun 3-Jul 2, 2002  '; 'Jul 3-30, 2002     '];
% Short version (for publications)
monthstr2 = ['Jun 2001'; 'Jul 2001'; 'Aug 2001'; 'Oct 2001'; ...
             'Nov 2001'; 'May 2002'; 'Jun 2002'; 'Jul 2002'];
monthstr = monthstr2;
% Compressed version (for filenames)
endstr  = ['Jun01'; 'Jul01'; 'Aug01'; 'Oct01'; ...
           'Nov01'; 'May02'; 'Jun02'; 'Jul02'];

GGall  = []; 
rhsall = [];
allsources = [];

for mindex = mlist,
  
  nr1  = depdata.Nr(mindex);   % nr1 is ignored
  ri   = rlist;
  rizn = ri;
  nrizn= length(rizn);
  risu = ri; nrisu = length(risu);
  risr = ri; nrisr = length(risr);
  
  % In fall the stability class is "slightly stable" (E), while
  % in summer we use "neutral" (D).
  if abs(mindex) == 4 || abs(mindex) == 5, stabclass = 'E'; end;
  stabclass = 'D';   % override: always use same stability class
  
  % These are dummy values used to calculate deposition fluxes only
  % (should be non-zero).
  source.Q0 = source.Qzn * 0 + 1;
  
  % Calculate the average depositions over all months
  % (for plotting purposes only).
  ritt = [1,4:8];    % leave out 2nd/3rd data sets (they're way off)
  meanDzn = mean(depdata.Zn(ritt,rizn),1);
  
  % Take the receptor data from the data file.
  if mindex > 0,
    recept.Dzn = depdata.Zn(mindex,rizn);
    recept.Dsr = depdata.Sr(mindex,risr);
    recept.Dsu = depdata.Su(mindex,risu);
  else
    mindex = abs(mindex);
    recept.Dzn = meanDzn;
    recept.Dsr = mean(depdata.Sr(ritt,risr),1);
    recept.Dsu = mean(depdata.Su(ritt,risu),1);
  end
  
  % Alternately, we can use the depositions (in kg) computed using
  % the "forward" code with source rates Q = [35, 80, 5, 5] T/yr.

  % For umin = 0.1, tskip = 20:
  %%recept.Dzn = 1.0e-04 * ...
  %% [ 0.177378268200547   0.708670618631528   0.191244296831422 ...
  %% 0.064945388392956   0.310790595942513   0.110541818040893 ...
  %% 0.040462421868453   0.056997414381165   0.107286103244816 ...
  %% 0.018884151161105   0.038896764809966   0.177378268200547 ...
  %% 0.310790595942513 ];

  % For umin = 0.1, tskip = 20:
  %% recept.Dzn = 1e-6 * ...
  %% [ 13.6559, 33.4813, 18.0888, 1.5989, 26.4213,  6.5377, 2.7183 ...
  %% 5.6984,  8.6759,  1.8340, 2.4882, 13.6559, 26.4213 ];
  %% recept.Dzn = recept.Dzn(rizn);

  % Add noise to the data (if necessary).
  recept.Dzn = recept.Dzn .* (1 + dnoise*(1-2*rand([1,nrizn])));
  recept.Dsu = recept.Dsu .* (1 + dnoise*(1-2*rand([1,nrisu])));
  recept.Dsr = recept.Dsr .* (1 + dnoise*(1-2*rand([1,nrisr])));
  
  % Generate a bar plot of the measured depositions.
  figure(1)
  tempdzn = zeros(1,nr);
  tempdsu = tempdzn; tempdsr = tempdzn;
  tempdzn(rizn) = recept.Dzn*1e6;  % in mg
  tempdsu(risu) = recept.Dsu*1e6;  % in mg
  tempdsr(risr) = recept.Dsr*1e9;  % in mu-g
  if nContam == 1,
    bar([tempdzn]', 'grouped')
  else                          
    bar([tempdzn;tempdsu;tempdsr]', 'grouped')
  end 
  if IsPlotMeanD, 
    hold on
    plot( rizn-0.22, meanDzn*1e6, 'ko' )
    hold off
  end
  set(gca,'XLim', [0.5,nr+0.5]);
  set(gca,'XTick',[1:nr]);
  set(gca, 'XTickLabel', recept.label)
  title( monthstr(mindex,:) );
  xlabel('Receptor'), ylabel('Amount deposited (mg)')
  if nContam >= 2,
    if IsPlotMeanD, 
      legend('Zn', 'SO_4', 'Sr*10^3', 'Zn Average'); 
    else
      legend('Zn', 'SO_4', 'Sr*10^3'); 
    end
  end
  grid on
  set(gca,'XGrid','off')
  if IsBWplots, 
    colormap(mybwcolmap)
    print('-depsc', ['invdep',endstr(mindex,:),'_bw.eps'])
  end
  colormap(mycolormap)
  set(gca,'YLim',[0,100])
  print('-depsc', ['invdep',endstr(mindex,:),'.eps'])
  shg 
  if IsPause, disp('Enter to continue ...'), pause, end 
  
  fprintf( 1, '\nINVERSE RUN %d:\n', mindex );
  disp( monthstr(mindex,:) );
  
  % Read the wind data for the appropriate month.
  switch mindex, 
   case {1,2,3}, wind = readwind( 'WindDataPartIIa.xls', 600, 0, umin ); % Jun'01 (for 1)
   case {4,5},   wind = readwind( 'WindDataPartIIb.xls', 600, 0, umin ); % Nov'01 (for 5)
   case 6,       wind = readwind( 'WindDataPartIIc.xls', 600, 0, umin ); % May'02
   case {7,8},   wind = readwind( 'WindDataJun3toJul2_2002B.xls', 600, 0, umin ); % Jun'02 (for 7)
   otherwise, 
    wind = readwind( 'WindDataJun3toJul2_2002B.xls', 600, 0, umin );
  end
  
  % Generate a radial bar histogram of the wind direction.  
  figure(2)
  windshift = pi; % Wind data has TO-angle, angular plots show FROM-angle
  wind2 = wind;   % for plotting purposes only
  wind2.dir = mod(wind2.dir+windshift, 2*pi);
  rose( wind2.dir )
  %%>> mmpolar( 'RTickAngle', 60, 'TTickValue', [0, 90, 180, 270], ...
  %%>>          'TTickLabel', {'E', 'N', 'W', 'S'} ); 
  %%>> %% set(gca, 'YTickLabelMode', 'manual')
  %%>> %% set(gca, 'YTickLabel', [])
  %%>> print('-depsc', ['invwind',endstr(mindex,:),'.eps'])
  %%>> shg
  %%>> if IsPause, disp('Enter to continue ...'), pause, end 
  
  plotwindrose( wind2, umin )
  print('-depsc', ['invwind2',endstr(mindex,:),'.eps'])
  shg
  if IsPause, disp('Enter to continue ...'), pause, end 
  
  % Now, solve the inverse problem.
  dt  = (wind.time(2) - wind.time(1)) * tskip;
  
  fprintf( 1, 'Constructing matrix (tskip=%d) ...\n', tskip );
  
  Gzn = zeros( nrizn, source.n );
  Gsr = zeros( nrisr, source.n );
  Gsu = zeros( nrisu, source.n ); % not used
  nwind = length(wind.dir);
  
  for k = [1 : tskip : nwind],
    
    k2 = min( k+tskip-1, nwind );
    if IsAvgWind,
      thetaavg = mean( wind.dir(k:k2) );  %% Doesn't average properly!!!
      Vwavg    = mean( wind.vel(k:k2) );
    else
      thetaavg = wind.dir(k);
      Vwavg    = wind.vel(k);
    end
      
    if mod(k,500) == 0, fprintf( 1, ' %5d ..', k ), end
    
    theta = +(pi/2 - thetaavg);
    Vw    = Vwavg;
    for i = 1 : source.n,
      xxx = (recept.x - source.x(i))*cos(theta) + (recept.y - source.y(i))*sin(theta);
      yyy =-(recept.x - source.x(i))*sin(theta) + (recept.y - source.y(i))*cos(theta);
      
      % Calculate Ermak solution with Q = 1, and accumulate into the G
      % matrices.  The units of G are s/m^3.
      warning( 'OFF', 'MATLAB:divideByZero' );
      Gzn(:,i) = Gzn(:,i) + ermak( xxx(rizn), yyy(rizn), recept.z(rizn), source.z(i), ...
                                   source.Q0(i), Vw, Vszns, Vdzns, stabclass )' / source.Q0(i);
      if nContam >= 2 && ~isempty(Gsr),
        Gsr(:,i) = Gsr(:,i) + ermak( xxx(risr), yyy(risr), recept.z(risr), source.z(i), ...
                                     source.Q0(i), Vw, Vssrs, Vdsrs, stabclass )' / source.Q0(i);
      end
      warning( 'ON', 'MATLAB:divideByZero' );
    end
    
  end
  fprintf( 1, '\n' );
  
  % Set up the linear least squares solve:
  %   * G is a (Nr x Ns)-matrix
  %   * rhs is a (Nr x 1)-vector
  %   * unknown Q is a (Ns x 1)-vector
  % where the number of sources (Ns) is less than 
  % the number of receptors (Nr).
  
  if nContam == 1,
    GG = Gzn;
    rhs = recept.Dzn' / (A * dt * Mzn * Vdzns); 
    % Set up positivity constraints
    Apos = -eye(4);
    bpos = [ 0; 0; 0; 0 ]; 
    % Set up the equality constraints
    Aeq = [ 0, 0, 1, -1 ]; 
    beq = [ 0 ];
  
  else

    if nContam == 2,
      GG  = [ Gzn, 0*Gzn; 0*Gsr, Gsr ];
      rhs = [ recept.Dzn/(Mzn*Vdzns), recept.Dsr/(Msr*Vdsrs) ]' / (A * dt); 
    else
      GG  = [ Gzn, 0*Gzn; 0*Gsr, Gsr; (Mso4/Mzn)*Gzn, (Mso4/Msr)*Gsr ];
      rhs = [ recept.Dzn/(Mzn*Vdzns), recept.Dsr/(Msr*Vdsrs), ...
              recept.Dsu/(Mso4*0.5*(Vdzns+Vdsrs)) ]' / (A * dt);  
    end
    
    % Set up negativity constraints as A * Q <= 0.
    Apos = -eye(8);
    npos = size(Apos,1);
    bpos = zeros(npos,1);
    
    % Set up the equality constraints.
    Aeq = [ 
    %%     Zinc         Strontium
        0, 0, 1,-1,    0, 0,  0,  0;
        0, 0, 0, 0,    0, 0,  1, -1;
        0, 0, 0, 0,    1, 0,  0,  0; 
        0, 0, 0, 0,    0, 1,  0,  0 ];
    if IsUseSr, Aeq = [ Aeq; 0, 0, 1, 0, 0, 0,-6e3, 0 ]; end
    %%if IsUseSr, Aeq = [ Aeq; 0, 0, 1, 0, 0, 0,-6e3, 0;
    %%                         0, 0, 0, 1, 0, 0, 0,-6e3 ]; end % dependent
    neq = size(Aeq,1);
    beq = zeros(neq,1);
  end
  
  if IsAllMonths,
    GGall = [GGall; GG];
    rhsall= [rhsall; rhs];
  end
  
  MM = [GG; Aeq];
  [UU,SS,VV] = svd(MM);
  fprintf( 1, 'Singular values:  [ ' );
  fprintf( 1, '%f  ', diag(SS) );
  fprintf( 1, ']\nCondition number: %f\n', cond(MM) );
  
  % Solve the linear system G*X=rhs using least squares.  
  % The solution X has units of mol/s.
  [X,res] = lsqlin( GG, rhs, Apos, bpos, Aeq, beq );
  %[X,res] = lsqlin( Gzn, rhs, Apos, bpos );
  %[X,res] = lsqlin( Gzn, rhs );
  
  % Fix up spurious non-zero values.
  X = max(X,0);
  
  % The solution Y has units of T/yr.
  fprintf(1, 'Zinc emission rates (T/yr):\n' );
  if nContam == 1,
    Y  = Mzn * X' / tpy2kgps;
    YY = Y;
    allsources(mindex).Qzn = YY;
    disp(Y);
    fprintf( '    Total: %e\n', sum(Y) );
  else
    allsources(mindex).Qzn = Mzn * X(1:end/2)'     / tpy2kgps;
    allsources(mindex).Qsr = Msr * X(end/2+1:end)' / tpy2kgps;
    allsources(mindex).Qsu = (Mso4/Mzn)*allsources(mindex).Qzn ...
                           + (Mso4/Msr)*allsources(mindex).Qsr;
    Y  = [ allsources(mindex).Qzn, ...
           allsources(mindex).Qsu, ...
           allsources(mindex).Qsr ];
    YY = reshape(Y, [4,3])';
    disp( YY(1,:) );
    fprintf( '    Total: %e\n', sum(YY(1,:)) );
    fprintf(1, 'Sulphate emission rates (T/yr):\n' );
    disp( YY(2,:) );
    fprintf( '    Total: %e\n', sum(YY(2,:)) );
    fprintf(1, 'Strontium emission rates (kg/yr):\n' );
    disp( YY(3,:)*1e3 );
    fprintf( '    Total: %e\n', sum(YY(3,:)) );
  end
  fprintf( 1, 'Residual = %e\n\n', res );
  
  % Plot the emission rate estimates as a bar plot.
  figure(3)
  if nContam >= 2, YY(3,:) = YY(3,:)*1e4; end
  bar(YY', 'grouped');
  if nContam >= 2, YY(3,:) = YY(3,:)/1e4; end
  set(gca,'XLim',[0.5,4.5]); 
  set(gca,'XTick',[1:4]);
  set(gca,'XTickLabel',['S1'; 'S2'; 'S3'; 'S4']);
  xlabel('Source'), ylabel('Emission rate (T/yr)')
  grid on
  set(gca,'XGrid','off')
  shg
  if IsPlotZnEst,
    hold on, plot( [1,2,3,4]-0.22, [35,80,5,5], 'ko' ), hold off
    legend('Zn    ','SO_4  ','Sr*10^4','Zn estimate' )
  else
    legend('Zn    ','SO_4  ','Sr*10^4')
  end
  title( monthstr(mindex,:) );
  if IsBWplots, 
    colormap(mybwcolmap)
    print('-depsc', ['invqrate',endstr(mindex,:),'_bw.eps'])
  end
  colormap(mycolormap)
  print('-depsc', ['invqrate',endstr(mindex,:),'.eps'])
  
  % Plot the total emission rates across all sources.
  figure(4)
  ysum  = sum(YY,2);
  ysum2 = ysum .* [1,1,1e4]';
  ymax  = max(ysum2(:));
  bar(ysum2)
  for i = 1 : length(ysum),
    text( i-0.18, ymax*0.04, num2str(ysum(i), 4) )
  end
  set(gca,'XTick',[1:3]);
  set(gca,'XTickLabel',['  Zn   '; ' SO_4  '; 'Sr*10^4']);
  title( monthstr(mindex,:) )
  ylabel('Total emission rate (T/yr)')
  grid on
  set(gca,'XGrid','off')
  shg
  if IsBWplots, 
    colormap(mybwcolmap)
    print('-depsc', ['invqtot',endstr(mindex,:),'_bw.eps'])
  end
  colormap(mycolormap)
  print('-depsc', ['invqtot',endstr(mindex,:),'.eps'])
  
end


if IsAllMonths,
  
  MM = [GGall;Aeq];
  [UU,SS,VV] = svd(MM);
  fprintf( 1, 'Singular values:  [ ' );
  fprintf( 1, '%f  ', diag(SS) );
  fprintf( 1, ']\nCondition number: %f\n', cond(MM) );
  
  fprintf( 1, '\nOPTIMIZE OVER ALL MONTHS:\n' );
  [X,res] = lsqlin( GGall, rhsall, Apos, bpos, Aeq, beq );
  X = max(X,0);
  
  fprintf(1, 'Zinc emission rates (T/yr):\n' );
  if nContam == 1,
    Y  = Mzn * X' / tpy2kgps;
    YY = Y;
    disp(Y);
    fprintf( '    Total: %e\n', sum(Y) );
  else
    Y = [ Mzn * X(1:end/2)', Mso4 * (X(1:end/2)' + X(end/2+1:end)'), ... 
          Msr * X(end/2+1:end)' ] / tpy2kgps;
    YY = reshape(Y, [4,3])';
    disp( YY(1,:) );
    fprintf( '    Total: %e\n', sum(YY(1,:)) );
    fprintf(1, 'Sulphate emission rates (T/yr):\n' );
    disp( YY(2,:) );
    fprintf( '    Total: %e\n', sum(YY(2,:)) );
    fprintf(1, 'Strontium emission rates (kg/yr):\n' );
    disp( YY(3,:)*1e3 );
    fprintf( '    Total: %e\n', sum(YY(3,:)) );
  end
  
  % Plot the emission rate estimates as a bar plot.
  figure(3)
  if nContam == 3, YY(3,:) = YY(3,:)*1e4; end
  bar( YY', 'grouped' );
  if nContam == 3, YY(3,:) = YY(3,:)/1e4; end
  set(gca,'XTick',[1:4]);
  set(gca,'XTickLabel',['S1'; 'S2'; 'S3'; 'S4']);
  xlabel('Source'), ylabel('Emission rate (T/yr)')
  grid on
  set(gca,'XGrid','off')
  shg
  if IsPlotZnEst,
    hold on, plot( [1,2,3,4]-0.22, [35,80,5,5], 'ko' ), hold off
    legend('Zn    ','SO_4  ','Sr*10^4','Zn estimate' )
  else
    legend('Zn    ','SO_4  ','Sr*10^4');
  end
  title( 'All months' )
  if IsBWplots, 
    colormap(mybwcolmap)
    print('-depsc', 'invqrateall_bw.eps')
  end
  colormap(mycolormap)
  print('-depsc', 'invqrateall.eps')
  
end

save 'invsave.mat'
toc

%% Ermak solution forward test
clear all
clc

profile on
tic

setparams;
setSolver;

[dep, dep2] = forwardErmak(source, recept, wind, xmesh,...
    ymesh, zmesh, dx, dy, dz, dt, tskip, nwind, Vszns, Vdzns, Mzn, stabclass,A,...
    IsAvgWind, IsOverDomain);

%% 
sourceUnit = source;
sourceUnit.Qzn = sourceUnit.Qzn./sourceUnit.Qzn;

%[depUnit, depUnit2] = forwardErmak(sourceUnit, recept, wind, xmesh,...
 %   ymesh, zmesh, dx, dy, dz, dt, tskip, nwind, Vszns, Vdzns, Mzn, stabclass,A,...
  %  IsAvgWind, IsOverDomain);

%depUnit

toc
profile viewer
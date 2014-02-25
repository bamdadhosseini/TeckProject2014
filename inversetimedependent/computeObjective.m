%% computeObjectiveAndGradient
%
% This function packages the forward solves to 
% return the value of the Tikhonov functional.
%

function [J,C, sensor] = computeObjective( sensor, source, wind,...
                    Vs, Vd, tskip, stabclass, R, alpha, targetQ)

    %% solve forward problem
    %
    % get accumulated concentrations
    C = getConcentrationAtReceptors(sensor, source, wind, Vs, Vd, stabclass, tskip);
    sensor = getAccumulatedConcentration( C, sensor, 1:length(sensor));
    %
    % get misfit
    sensor = getMeasurementAtSensor( sensor );
    sensor = getMisfitAtSensor( sensor );
    %
    % evaluate Tikhonov Functional
    J = evalTikhonovFunctional( source, sensor, R, alpha, targetQ);
    
end
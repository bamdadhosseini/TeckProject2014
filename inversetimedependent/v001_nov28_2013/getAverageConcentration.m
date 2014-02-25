%% getAverageConcentration( sensor, M, dt )


function sensor = getAverageConcentration( C, sensor, M, dt )

% loop over sensors 

for iSens = 1:length(sensor)

    thissens  = sensor(iSens); % helps with a cleaner code, not a big price to pay
    thisconc = C(iSens, :);   % again just to make things more clear
    
    % partition the concentration (mol/m^s) into a matrix where each column
    % is the concentration during measurement period 
    
    data = [];
    for iMeas = 1:length(thissens.start_indx)
       
        data(:,iMeas) = thisconc(thissens.start_indx(iMeas):thissens.end_indx(iMeas));
       
    end
    %size(C)
    %size(data)
    
    % now compute the averages
    
    % number of averaging steps between each start and end time
    numberofaverages = ceil((thissens.end_indx(1)- thissens.start_indx(1))/thissens.ave_period);
    %disp(['numberofaverages ', num2str(numberofaverages)])
    aveConc = zeros(numberofaverages,size(data,2));
    
    %disp(numberofaverages)
    
    for iAve = 1:numberofaverages
     %   disp(['iAve ', num2str(iAve)]);
     %   disp(['startindx: ', num2str(1+(iAve-1)*thissens.ave_period), ' to ', ...
     %       num2str(iAve*thissens.ave_period) ] );
        sensor(iSens).forwardmeasurement(iAve,:) = (1/thissens.ave_period)*...
                        sum(data(1+(iAve-1)*thissens.ave_period:iAve*thissens.ave_period,:), 1);
     %   pause;
    end
end

end
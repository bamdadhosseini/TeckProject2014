%% getAverageConcentation( sensor, M, dt )


function sensor = getAccumulatedConcentration( C, sensor, index )

% loop over sensors 

for iSens = index

    thissens  = sensor(iSens); % helps with a cleaner code, not a big price to pay
    thisconc = C(iSens, :);   % again just to make things more clear
    
    % partition the concentration (mol/m^s) into a matrix where each column
    % is the concentration during measurement period 
    
    %iSens
    %thissens
    %pause
    %thisconc(thissens.start_indx(1):thissens.end_indx(1))
    %pause
    %disp(iSens)
    %size(thisconc)
    data = [];
    for iMeas = 1:length(thissens.start_indx)
        %disp(iMeas)
        %disp(thissens.start_indx(iMeas))
        %disp(thissens.end_indx(iMeas))
     
        data(:,iMeas) = thisconc(thissens.start_indx(iMeas):thissens.end_indx(iMeas));
       
        %disp(size(data(:,iMeas)));
        %pause
    end

    % now compute the averages
    
    % number of averaging steps between each start and end time
    numberofaverages = max(1,ceil((thissens.end_indx(1)- thissens.start_indx(1)+1)/thissens.accum_period));
    %disp(['numberofaverages ', num2str(numberofaverages)])
    aveConc = zeros(numberofaverages,size(data,2));
    
    %disp(numberofaverages)
    %pause
    for iAve = 1:numberofaverages
        %disp(['iAve  ', num2str(iAve)])
        %pause
        sensor(iSens).accum_concentration(iAve,:) = ...
                        sum(data(1+(iAve-1)*thissens.accum_period:iAve*thissens.accum_period,:), 1);
    end
end

end
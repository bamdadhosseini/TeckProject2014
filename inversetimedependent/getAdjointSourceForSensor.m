%% getAdjointSourceForSensor
%
% construct the source term of the Adjoint equation using the 
% computed misfit and measurement data such as start and 
% end times as well as averaging periods.

function AdjSource = getAdjointSourceForSensor( sensor, sizeAdjC )


    % initialize the matrix of the source term, this is 
    % a matrix of size R X N where R is the number of sensors 
    % and N is the number of time steps (same size as the 
    % forward concentration. Each column is the source values 
    % at the corresponding time step. So AdjC(i,j)
    % is the source term for sensor(i) at time t_j
    
    AdjSource = zeros(sizeAdjC);

    %disp(['size Adjoint source :', num2str(size(AdjC))]);
    % iterate over each sensor, each one will have a 
    % different source term 
    
    for iSens = 1:length(sensor)
       
        % the source term of the adjoint equation for each 
        % sensor is constant during each averaging period 
        % and equal to the corresponding misfit and zero 
        % outside of measurement times.
        
        % loop over start times of measurements
        
        thissens = sensor(iSens); % to keep the code more readable
        %disp(['sensor : ', num2str(iSens)])
        %disp(['sizeAdjC : ', num2str(sizeAdjC)])
        for k = 1:size(thissens.start_indx,2)
           %disp(['k : ',num2str(k)]);
           
           % number of averaging steps between each start and end time
           %numberofaverages = ceil((thissens.end_indx(1)- thissens.start_indx(1))/thissens.accum_period);
            
           % repeat the misfits for the same number of times as 
           % the averaging periods 
           
           MM = repmat(thissens.misfit(:,k).*thissens.scale.*thissens.weight,...
               1, thissens.accum_period)';
           %MM = MMM(:)';
          % MM(end:sizeAdjC(2)) = 0; % extend the end by zero 
                                    % if we didnt have enough time to take 
          %                          % measurements.
           
           %disp(['size MMM:', num2str(size(MMM))]);
%            disp(['size MM(:):', num2str(size(MM(:)'))]);
%            disp(['startandendtime: ', num2str([thissens.start_indx(k),  thissens.end_indx(k)])]);
%            disp(size(AdjSource(iSens,thissens.start_indx(k):thissens.end_indx(k))))
%            pause
           AdjSource(iSens,thissens.start_indx(k):thissens.end_indx(k)) = MM(:)';
        end 
    end
end
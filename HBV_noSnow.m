function Q = HBV_noSnow(Params_struct, ClimateData, doMonthly) 
%HBV_noSnow runs the HBV rainfall runoff model without any routing.
%   
%   Adapted from 
%
%   See AghaKouchak A., Habib E., 2010, Application of a Conceptual 
%   Hydrologic Model in Teaching Hydrologic Processes, International 
%   Journal of Engineering Education, 26(4), 963-973.
% 
%   AghaKouchak A., Nakhjiri N., and Habib E., 2013, An educational model 
%   for ensemble streamflow simulation and uncertainty analysis, Hydrology 
%   and Earth System Sciences, 17, 445-452, doi:10.5194/hess-17-445-2013.
%
%   Written by Andrew Western, University of Melbourne, March 2018


   
    % Extract climate data.
    precip = ClimateData(:,4);
    PET = ClimateData(:,5);

    % Extract parameters
    fc = Params_struct.fc; % maximum soil water storage
    beta = Params_struct.beta; % runoff non-linearity
	pwp = fc * Params_struct.pwp; % soil water AET control
    l = Params_struct.l; % fast response threshold
    k0 = Params_struct.k0; % fast response recession
    k1 = Params_struct.k1; % store 1 recession
    kp = Params_struct.kp; % store 1 percolation recession
    k2 = Params_struct.k2; % store 2 recession
    

    % Set initial conditions
    soil=0.3*pwp;      % Initial soil moisture store
    s1 = 0;% Initial routing store
    s2 = 0; % Initial groundwater store

    Q=nan(length(precip),1);

    % Start DAILY LOOP
    for z=1:length(precip);
        % Actual ET
        if(soil>pwp)
            aet=min(PET(z),soil);
        elseif soil > 0
            aet=min(PET(z)*(soil/pwp),soil);
        else
            aet = 0;
        end
        
        % Runoff and soil storage
        directRunoff = precip(z)*(soil/fc)^beta;
        soil = soil+precip(z)-directRunoff-aet;
        if soil > fc
            directRunoff = directRunoff + soil - fc;
            soil = fc;
        end
        
        % Stores
        percolation = s1*kp;
        recession_0 = max(0,s1-l)*k0;
        recession_1 = s1*k1;
        recession_2 = s2*k2;
        
        s1 = s1 + directRunoff - recession_0 - recession_1 - percolation;
        s2 = s2 + percolation - recession_2;
        Q(z)=recession_0 + recession_1 + recession_2;
            
    end
    % End DAILY LOOP
    
    % Add year, month and day to flow vector.
    Q = [ ClimateData(:,1:3), Q];
    
    % If the calibration is to be undertaken using monthly flow data,
    % aggregate the flow dat to monthly.
    if doMonthly        
        Q = convertDailyToMonthly(Q);       
    end
end

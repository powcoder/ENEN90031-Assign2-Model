function [Params_best, bestObjFn] = ...
    SCE_calibration(ParamNames2Calib, ParamsInitial, Params_LowerBound, ...
    Params_UpperBound, timeStart, timeEnd, climateData, observedFlow, ...
    doMonthly, optimisationSettings)
% SCE_CALIBRATION undertakes a shuffled complex evolution calibration of the 
% HBV rainfall runoff model. NOTE: the SCE setting are and the parameters
% to calibrate are hard-coded into the function (see lines 34-51). The SCE 
% settings can be edited.
%
% Input summary
%
% ParamsInitial:a structure with parameter values. Importantly, if a parameter is
%               not to be calibrated then the parameter value within this structure
%               will be used.			
%
% Params_LowerBound: a structure with parameter value at the lower boundaries.
%
% Params_UpperBound: a structure with parameter value at the upper boundaries.
%
% timeStart and timeEnd - the start and end time of the simulation in
%               matlab datenumber format.
%
% ClimateData: 	daily climate data in the following format: year, month,
%               day, precip (mm/day), potential evap (mm/day).
%
% observedFlow:daily observed streamflow in the following format: year, month,
%              	day, flow (mm/day).
%
% doMonthly: 	logical scaler for turning on calibration to monthly
%            	observation data.
%
% Output summary
% Params_best:  a structure with parameter values at the lowest objective 
%               function value.
% bestObjFn:    the lowest objective function value.
%
%
% Example:
%   timeStart = datenum(sevensCreekClimate(1,1), sevensCreekClimate(1,2), sevensCreekClimate(1,3));
%   timeEnd = datenum(sevensCreekClimate(end,1), sevensCreekClimate(end,2), sevensCreekClimate(end,3));
%   doMonthly = false;
%
%   [Params_best, bestObjFn] = ...
%     SCE_calibration( Params_initialModel, Params_LowerBound, ...
%     Params_UpperBound, timeStart, timeEnd, sevensCreekClimate, sevensCreek_mm_day, doMonthly);
%
% Written by Tim Peterson, University of Melbourne May, 2010 and updated in
% May 2013. Example added in 2016.
    
    % Set the names of the parameters to calibrate.
    % If a parameter is not to be calibrated  
    % then it should be DELETED.    
%     ParamNames2Calib = {'fc'; 'beta'; 'pwp'; 'l'; 'k0'; 'k1'; 'kp'; 'k2'};
               
    % Set calibration options.
    % Please modify these to invesigate if SCE converges to the global
    % minimum!
    maxn = optimisationSettings.maxn;
    kstop = optimisationSettings.kstop;    
    pcento = optimisationSettings.pcento;    
    peps = optimisationSettings.peps;
    ngs = optimisationSettings.ngs;
    iseed = floor(rand(1)*100000);
    iniflg =  0;
    
    % If the calibration is to be undertaken using monthly flow data,
    % aggregate the flow dat to monthly.
    
    % Extract initial parameter values. 
    for i=1: length(ParamNames2Calib)
        params(i) = ParamsInitial.(ParamNames2Calib{i});
    end
    
    % Extract upper and lower bounds parameter values. 
    for i=1: length(ParamNames2Calib)
        params_lower(i) = Params_LowerBound.(ParamNames2Calib{i});
    end    
    for i=1: length(ParamNames2Calib)
        params_upper(i) = Params_UpperBound.(ParamNames2Calib{i});
    end    
    
    % Check that the lower bound parameter values are less than the upper
    % bound values.
    if any( params_upper - params_lower < 0)
        error('Lower bound parameter values cannot be greater than the upper bound values!');
    end
    
    % Do SCE calibration using the objective function SSE.m 
    [bestParams_temp, bestObjFn] = sceua(@objectiveFunction, params, params_lower, params_upper, maxn, kstop, pcento, peps, ngs, iseed, iniflg, ...
    ParamNames2Calib, ParamsInitial, timeStart, timeEnd, climateData, observedFlow, doMonthly);    
    
    % Assign the initial parameter values back to the parameters data structure.
    Params_best = ParamsInitial;   
    for i=1: length(ParamNames2Calib)
        Params_best.(ParamNames2Calib{i}) = bestParams_temp(i);
    end

end
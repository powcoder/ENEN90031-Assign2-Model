function [obj, residuals, Q] = objectiveFunction( parameterVector, parameterNames,  ...
    modelParameters, timeStart, timeEnd, climateData, observedFlow, doMonthly)
%OBJECTIVEFUNCTION Run the model and calculate an objective function.
%
% This sets the model parameters to those supplied, runs the model for the
% specified period and calculates an objective function value plus model
% residuals.
%
% Syntax:
%   [obj, residuals] = objectiveFunction( parameterVector, parameterNames,  ...
%       modelParameters, timeStart, timeEnd, climateData, observedFlow, doMonthly)
%
% Inputs:
%   parameterVector - nx1 vector of parameter values for to be used in the
%   model run.
%
%   parameterNames - cell vector containing the name of the parameters to
%   set. The parameter names must be listed within the inputs
%   'parameterStructure'.
%   
%   modelParameters - Structural variable containing the parameters for
%   the model. 
%
%   timeStart and timeEnd - the start and end time of the simulation in
%   matlab datenumber format
%
%   climateData - a matrix with row for each day containing 
%   year, month, day, rainfall, and potential evapotranspiration values
%
%   observedFlow - a matrix with row for each day (or month) containing 
%   year, month, day, observed flow values
%
%   doMonthly - if set to true monthly flows are used for the objective
%   function, residuals and discharge.  If true, observedFlow must be
%   monthly with the day for each month set to the last day of the month
%
%   
%
% Outputs:
%   obj - a scalar containing the objective function value.
%
%   residuals - an nx1 vector of model residual values (simulated -
%   observed)
%
%   Q - simulated flows in mm for each day (or month if doMonthly = true)
%
% Example:
%   
%   % Set the parameters to be changed
%   parameterNames = {'INSC'; 'SMSC'};
% 
%   % Set the values for those parameters
%   parameterVector = [2, 150];
% 
%   % Set times
%   timeStart = datenum(BrokenRiver_Climate_Daily_mm(1,1), ...
%       BrokenRiver_Climate_Daily_mm(1,2), BrokenRiver_Climate_Daily_mm(1,3));
%   timeEnd = datenum(BrokenRiver_Climate_Daily_mm(end,1), ...
%       BrokenRiver_Climate_Daily_mm(end,2), BrokenRiver_Climate_Daily_mm(end,3));
% 
%   % Run the model and obtain the objective function value.
%   [obj, residuals] = objectiveFunction( parameterVector, parameterNames,  ...
%       Params_initialModel, timeStart, timeEnd, BrokenRiver_Climate_Daily_mm, ...
%       BrokenRiver_Flow_Daily_mm, false);
%

modelParameters = setParameterValues( parameterVector, parameterNames, modelParameters);
spinUp = 1; %in years

% Trim climate data
theDates = datenum(climateData(:,1), climateData(:,2), climateData(:,3));
if(theDates(1) > timeStart)
    error('Climate data begins after start time');
    return
end
if(theDates(end) < timeEnd)
    error('Climate data finishes before end time');
    return
end

tstart = datevec(timeStart);
spinupStart = datenum(tstart(1)-spinUp,tstart(2),tstart(3));
numSpinupSteps = timeStart - spinupStart;

if theDates(1) <= spinupStart;
    % real forcing available for spin up
    filt = theDates >= spinupStart & theDates <= timeEnd;
    climateData = climateData(filt,:);
else
    % reuse early data for spinup
    filt = theDates >= timeStart & theDates <= timeEnd;
    climateData = climateData(filt,:);
    climateData = [climateData(1:numSpinupSteps,:); climateData];
    % replace spinup  dates
    temp = datevec(spinupStart:timeStart-1);
    climateData(1:numSpinupSteps,1:3) = temp(:,1:3); 
    warning('Insufficient spinup data, using first part of modelling period for spin up');
end

% Run the model
Q = HBV_noSnow(modelParameters, climateData, false);

% Remove spinup period
Q = Q(numSpinupSteps+1:end,:);

% convert observations to monthly if doMonthly true
if doMonthly
    Q= convertDailyToMonthly(Q);
    observedFlow = convertDailyToMonthly(observedFlow);
end

% get observation dates
obsDates = datenum(observedFlow(:,1), observedFlow(:,2), observedFlow(:,3));


% Calculate the residuals
% Trim series to match start/end dates before calculating residuals
modDates=datenum(Q(:,1),Q(:,2),Q(:,3));
firstDate = max([modDates(1),obsDates(1), timeStart]);
lastDate = min([modDates(end),obsDates(end), timeEnd]);
Qfilt = modDates >= firstDate & modDates <= lastDate;
obsFilt = obsDates >=firstDate & obsDates <= lastDate;

residuals = Q(Qfilt,4) - observedFlow(obsFilt,4);
residuals = residuals(~isnan(residuals));

% Calculate the objective function
obj = residuals' * residuals;

end
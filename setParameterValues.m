function parameterStructure = setParameterValues( parameterVector, parameterNames, parameterStructure)
%SETPARAMETERVALUES Input a vector of parameter values to a structure variable.
%
% This assigns an input vector of parameter values to a structure variable.
% The function is used to input parameter values from calibration inot a
% model structure.
%
% Syntax:
%   parameterStructure = setParameterValues( parameterVector, ...
%   parameterNames, parameterStructure)
%
% Inputs:
%   parameterVector - nx1 vector of parameter values for those being
%   calibrated.
%
%   parameterNames - cell vector containing the name of the parameters to
%   calibrate. The parameter names must be listed within the inputs
%   'parameterStructure'.
%   
%   parameterStructure - Structural variable containing the parameters for
%   the model. 
%
% Outputs:
%   parameterStructure - Updated structural variable containing the new 
%   parameter values input to the function.
%
% Example:
%   
%   % Set the parameters to calibrate
%   parameterNamesForCalib = {'hydroCapacity'; 'pumpCapacity'};
%
%   % Set parameter values to be input to the model structure.
%   parameterFromCalib = [ 3.2; 3.1];
%
%   % Update structure variable with  parameter values.
%   sysInfo_new = setParameterValues( parameterFromCalib, ...
%   parameterNamesForCalib, sysInfo)
%

    % Check the inputs are of the correct form.
    if ~isnumeric(parameterVector);
        error('The input objectiveFunction "parameterVector" must be a vector of parameter values.');
    elseif ~isvector(parameterVector) && ~isvector(parameterNames)
        error('The input objectiveFunction "parameterVector" and "parameterNames" must be vectors.');        
    elseif ~ischar(parameterNames) && ~iscell(parameterNames)
        error('The input objectiveFunction "parameterNames" must be a single parameter name or a cell vector of parameter names.');
    elseif length(parameterNames) ~=length(parameterVector)
        error('The input objectiveFunction "parameterVector" and "parameterNames" must be vectors of the same size.');
    end
    
    % Get the number of parameters to calibrate.
    numParams = length(parameterNames);
    
    % Cycle through each parameter name, get the value from the parameter
    % vector and assign it to the parameter structure variable.
    for i=1:numParams
        parameterStructure.(parameterNames{i}) = parameterVector(i);
    end
end


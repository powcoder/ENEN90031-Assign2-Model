function parameterVector = getParameterValues( parameterNames, parameterStructure)
%GETPARAMETERVALUES Get vector of parameter values from a structure variable.
%
% This retrives a vector of parameter values from a structure variable.
%
% Syntax:
%   parameterVector = getParameterValues( parameterNames, parameterStructure)
%
% Inputs:
%   parameterNames - cell vector containing the name of the parameters to
%   calibrate. The parameter names must be listed within the inputs
%   'parameterStructure'.
%   
%   parameterStructure - Structural variable containing the parameters for
%   the model. 
%
% Outputs:
%   parameterVector - nx1 vector of parameter values for those defines in 
%   'parameterNames'.
%
% Example:
%   
%   % Set the parameters to calibrate
%   parameterNamesForCalib = {'hydroCapacity'; 'pumpCapacity'};
%
%
%   % Update structure variable with  parameter values.
%   parameterValues = getParameterValues(parameterNamesForCalib, sysInfo)
%

    % Check the inputs are of the correct form.
    if ~isvector(parameterNames)
        error('The input objectiveFunction "parameterVector" and "parameterNames" must be vectors.');        
    elseif ~ischar(parameterNames) && ~iscell(parameterNames)
        error('The input objectiveFunction "parameterNames" must be a single parameter name or a cell vector of parameter names.');
    end
    
    % Get the number of parameters to calibrate.
    numParams = length(parameterNames);
    
    % Cycle through each parameter name, get the value from the parameter
    % vector and assign it to the parameter structure variable.
    for i=1:numParams
        parameterVector(i,1) = parameterStructure.(parameterNames{i});
    end
end


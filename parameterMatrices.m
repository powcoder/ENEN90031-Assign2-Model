function M = parameterMatrices(parameters, parametersLower, parametersUpper, numRealisations, parameterNames)
%parameterMatrices generates matrices of stochastic parameters
%   
%   This is for generating a matrix of parameter values given a lower and
%   upper bound. Stochastic parameters are generated for
%   parameterNames and values are taken from parameters for the rest.  
%   numRealisations parameter sets are generated from a uniform distribution.  
%
%   Written by Andrew Western, University of Melbourne, May 2014


%start with Uniform [0,1] distribution
% apply ranges to variables by linear scaling from [0,1] range to [min,max] range.
paramNames = fieldnames(parameters);
numParams=length(paramNames);

M=NaN(numRealisations,numParams);

for i=1:numParams
    M(:,i)=parameters.(paramNames{i});
end

for j=1:length(parameterNames)
    %find the correct column
    for i=1:numParams
        if strcmp(paramNames(i),parameterNames(j)) 
            break
        end
    end
    M(:,i)=unifrnd(parametersLower.(parameterNames{j}), parametersUpper.(parameterNames{j}), numRealisations,1);
end

end


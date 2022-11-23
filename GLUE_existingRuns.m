function [ medianPrediction, predictionCIs, parameterDistributions, behaviouralParams, ...
behaviouralPredictions, likelihood, pcntObsAboveUpperCI, pcntObsBelowLowerCI ] = ...
GLUE_existingRuns( ParamNames4Glue, parameterValues, performance, simulated, threshold, ...
observed, doMonthly, forceToZero, effectiveZero )

% GLUE_EXISTINGRUNS: This function undertakes a GLUE analysis of a 
%   rainfall-runoff model.
%
% [ medianPrediction, predictionCIs, parameterDistributions, behaviouralParams, ...
% behaviouralPredictions, likelihood, pcntObsAboveUpperCI, pcntObsBelowLowerCI ] = ...
% HBV_GLUE_existingRuns( ParamNames4Glue, parameterValues, sse, simulated, threshold, observed )
%
% Input summary
% ParamNames4Glue: a cell array with the parameter names 
% parameterValues: a matrix of parameter vlues - rows corresponde to realisations,
%               columns correspond to parameters.
% performance:          sum of squared error for each model run
% simulated:    a matrix of simulated flows, rows correspond to dates and
%               columns correspond to realisations
% threshold:    the GLUE threshold for a simulation to be accepted into the
%               uncertainty estimation. It is defined as the ratio of the
%               sum of square of errors divided by the observed sum of
%               squares.
% observed:     daily observed streamflow in the following format: year, month,
%              	day, flow.
% doMonthly: 	logical scalr for turning on calibration to monthly
%            	observation data.
% forceToZero:  if true,very small values of lower confidence bound are 
%               foreced to zero
% effectiveZero:threshold value below which lower confidence bound is
%               assumed to be zero.  Appropriate value would be lowest
%               practically measurable flow.
%
% Output summary:
% medianPrediction: a matrix with the first column equalling the time (in years) and the
%               second column equaling the median lake phosphorus concentration at the start of
%               that time point.
% predictionCIs:a matrix with time, 5% and 95% confidence limits in
%               columns.
% parameterDistributions: a structure of arrays with col 1 = parameter
%               value, col 2 = cumulative likelihood.
% behaviouralParams: an array (columns are s, r, m, q) of all behavioural
%               parameter sets.
% behaviouralPredictions: an array where each row is a timeseries of
%               predicted P for all behavioural runs
% likelihood:   a vector of likelihoods corresponding to behavioural runs
% pcntObsAboveUpperCI: fraction of observation above the calculated upper
%               confidence interval.
% pcntObsBelowLowerCI: fraction of observation below the calculated lower
%               confidence interval.
%
% Example
%
%  threshold = 0.5;
%
% [ medianPrediction, predictionCIs, parameterDistributions, behaviouralParams, ...
% behaviouralPredictions, likelihood, pcntObsAboveUpperCI, pcntObsBelowLowerCI ] = ...
% HBV_GLUE_existingRuns( parameterNames, paramsForSimulated, sseSimulated, ...
% simulated, threshold,  ovensFlow_mmperDay );
%  
%
%
% Written by Tim Peterson and Andrew Western, University of Melbourne May, 
% 2013.

    display(['Staring GLUE analysis ...']);
    
    display(['   Filtering the model runs by the threshold value of ',num2str(threshold),' ...']);
    
    % Find length of runs
    simulated = simulated';
    numberOutputTimes = length(simulated(1,:));
    simDates = datenum( observed(:,1), observed(:,2), observed(:,3));
    realisations = length(parameterValues(:,1));
    % Find behavioural runs
    observedDataFilter = ~isnan(observed(:,end));
    meanObs = mean(observed(observedDataFilter,end));
    observedSumSquares = (observed(observedDataFilter,end)'-meanObs)* ...
        (observed(observedDataFilter,end)-meanObs);
    behaviouralFilter=performance<=threshold;
    
    if sum(behaviouralFilter) <20
        error('There were <20 model simulations within the threshold. Try increasing the realistions or the threshold');
    end
    
    numberBehavioural = sum(behaviouralFilter);
    behaviouralNormalisedSSE = performance(behaviouralFilter);
    behaviouralParams = parameterValues(behaviouralFilter,:);
    behaviouralPredictions = simulated(behaviouralFilter,:);
        
    display('   Calculatting the likelihoods and parameter distributions ...');
    % Calculate the likelihoods
    likelihood = 1 - behaviouralNormalisedSSE;
    min_lh = min(min(likelihood),0);
    likelihood = likelihood - repmat(min_lh,numberBehavioural,1);
    sum_lh = sum(likelihood);
    likelihood = likelihood ./ repmat(sum_lh,numberBehavioural,1);
    
    % Remove the zero likelihood case
    filt = likelihood ~= 0;
    likelihood = likelihood(filt);
    behaviouralNormalisedSSE = behaviouralNormalisedSSE(filt);
    behaviouralParams = behaviouralParams(filt,:);
    behaviouralPredictions = behaviouralPredictions(filt,:);
    
    % Calculate parameter distributions
    for i=1: length(ParamNames4Glue)
        parameterDistributions.(ParamNames4Glue{i}) = sortrows([behaviouralParams(:,i),likelihood]);
        parameterDistributions.(ParamNames4Glue{i})(:,2) = cumsum(parameterDistributions.(ParamNames4Glue{i})(:,2));
    end
    
    
    display('   Calculating the simulation confidence intervals...');
    % calculate 90% confidence intervals for  predictions
    predictionCIs = zeros(numberOutputTimes,3);
    medianPrediction = zeros(numberOutputTimes,2);
    
    % get the output times
    predictionCIs(:,1) = simDates;
    medianPrediction(:,1) = simDates;
    for i = 1:numberOutputTimes
        temp = sortrows([behaviouralPredictions(:,i),likelihood]);
        temp(:,2) = cumsum(temp(:,2));
        medianPrediction(i,2) = interp1(temp(:,2),temp(:,1),0.5);
        predictionCIs(i,2) = interp1(temp(:,2),temp(:,1),0.05);
        predictionCIs(i,3) = interp1(temp(:,2),temp(:,1),0.95);
    end
    
    % force small to zero
    if forceToZero
                predictionCIs(predictionCIs(:,2)<effectiveZero,2) = 0;
    end

    % Transpose so time is by rows and realisations by columns.
    behaviouralPredictions = behaviouralPredictions';
    
    % Plot the simulated P
    figure()    
    plot(simDates, observed(:,end),simDates, medianPrediction(:,2),simDates, predictionCIs(:,2),simDates, predictionCIs(:,3)); 
    legend('Observed','Modelled-Median', 'Modelled-5th', 'Modelled-95th');
    title(['GLUE Prediction uncertainty for ',num2str(realisations),' realsiation and a threshold of ', num2str(threshold)]);
    xlabel('Time (years)');
    if doMonthly
        ylabel('Streamflow (mm/month)');        
    else
        ylabel('Streamflow (mm/d)');        
    end
    box on;    
    
    
    pcntObsBelowLowerCI = sum(observed(:,end) < predictionCIs(:,2))/numberOutputTimes;
    pcntObsAboveUpperCI = sum(observed(:,end) > predictionCIs(:,3))/numberOutputTimes;
    display('   Finished GLUE anaalysis.');
    display(['   The fraction of observations above the upper confidence interval =',num2str(pcntObsAboveUpperCI)]);
    display(['   The fraction of observations below the lower confidence interval =',num2str(pcntObsBelowLowerCI)]);
    
    
end

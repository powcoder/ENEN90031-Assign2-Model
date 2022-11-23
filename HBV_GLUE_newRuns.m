function [ medianPrediction, predictionCIs, parameterDistributions, behaviouralParams, ...
behaviouralPredictions, likelihood, pcntObsAboveUpperCI, pcntObsBelowLowerCI, ...
MCparameters, performance, simulated, outputDates] = ...
HBV_GLUE_newRuns( start_year, end_year, parameters, parametersLowerBound, parametersUpperBound, ...
parameterNames, ClimateData, observed, doMonthly, realisations, threshold, forceToZero, effectiveZero )

% This script undertakes a GLUE analysis of the 
% HBV rainfall runoff model.
%
% [ medianPrediction, predictionCIs, parameterDistributions, behaviouralParams, ...
% behaviouralPredictions, likelihood, pcntObsAboveUpperCI, pcntObsBelowLowerCI, ...
% MCparameters, performance, simulated, outputDates] = ...
% HBV_GLUE_newRuns( start_year, end_year, parameters, parametersLowerBound, parametersUpperBound, ...
% parameterNames, ClimateData, observed, doMonthly, realisations, threshold )
% 
% Input summary
% start_year:   simulation start year 
% end_year:     simulation start year 
% parameters:a structure with parameter values. Importantly, if a parameter is
%               not to be calibrated then the parameter value within this structure
%               will be used.			
% parametersLowerBound: a structure with parameter value at the lower boundaries.
% parametersUpperBound: a structure with parameter value at the upper boundaries.
% ClimateData: 	daily climate data in the following format: year, month,
%               day, precip (mm/day), potential evap (mm/day).
% Observed_Flow:daily observed streamflow in the following format: year, month,
%              	day, flow (mm/day).
% doMonthly: 	logical scalr for turning on calibration to monthly
%            	observation data.
% realisations: the number of random parameter sets as an integer.
% threshold:    the GLUE threshold for a simulation to be accepted into the
%               uncertainty estimation. It is defined as the ratio of the
%               sum of square of errors divided by the observed sum of
%               squares.
% forceToZero:  forces very small values of lower confidence bound to zero
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
% MCparameters: the parameters used for the full set of runs
%
% performance:  the sse of full set of runs
%
% simulated:    a matrix of simulated flows for runs where sse < observed
%               sum of squared deviation from the mean
%
% outputDates: 	the dates corresponding to output times
% 
% paramsForSimulated: the parameters for the runs in simulated
%
% sseSimulated: the sum of sqared error for the runs in simulated
%
%
% Example
%
%  start_year = 1981;
%  end_year = 2006;
%  doMonthly = false;
%  realisations = 10000;
%  parameterNames = fieldnames(Params_initialModel);
%  threshold = 0.5;
%  
%
% [ medianPrediction, predictionCIs, parameterDistributions, behaviouralParams, ...
%   behaviouralPredictions, likelihood, pcntObsAboveUpperCI, pcntObsBelowLowerCI, ...
%   MCparameters, performance, simulated, outputDates] = ...
%   HBV_GLUE_newRuns( start_year, end_year, Params_initialModel, Params_LowerBound, Params_UpperBound, ...
%   parameterNames, sevensCreekClimate, sevensCreek_mm_day, doMonthly, realisations, threshold );
% 
%
% Written by Tim Peterson and Andrew Western, University of Melbourne May, 
% 2013.

    display('Staring GLUE analysis ...');

    % Trim to min of start climate and start obs flow data.
    ClimateDates = datenum( ClimateData(:,1), ClimateData(:,2), ClimateData(:,3));
    FlowDates = datenum( observed(:,1), observed(:,2), observed(:,3));
    startDate = max( max( min(ClimateDates) , min(FlowDates)), datenum(start_year, 1, 1) );
    endDate = min( min( max(ClimateDates) , max(FlowDates)), datenum(end_year, 12, 31) );
    ClimateDates_filt = ClimateDates >= startDate & ClimateDates <= endDate;
    FlowDates_filt = FlowDates >= startDate & FlowDates <= endDate;    
        
    % set up parameter matrix
    MCparameters = parameterMatrices(parameters, parametersLowerBound, parametersUpperBound, realisations, parameterNames);

    % run the model and calculate performance measures for grouping runs
    [performance, simulated, outputDates ] = HBV_MonteCarlo(MCparameters, ...
        parameters, startDate, endDate, ClimateData, observed, doMonthly, true, true); 
    
    if doMonthly
        observed = convertDailyToMonthly(observed);
        FlowDates = datenum( observed(:,1), observed(:,2), observed(:,3));
        FlowDates_filt = FlowDates >= startDate & FlowDates <= endDate;    
    end
    observedFilt = observed(FlowDates_filt,:);        

    % Find length of runs
    simulated = simulated';
    numberOutputTimes = length(simulated(1,:));
    simDates = datenum( outputDates(:,1), outputDates(:,2), outputDates(:,3));
    
    % Find behavioural runs
    observedDataFilter = ~isnan(observed(:,end));
    behaviouralFilter=performance<=threshold;
    
    numberBehavioural = sum(behaviouralFilter);
    if numberBehavioural <20
        error('There were <20 model simulations within the threshold. Try increasing the realistions or the threshold');
    end    
    behaviouralNormalisedSSE = performance(behaviouralFilter);
    behaviouralParams = MCparameters(behaviouralFilter,:);
    behaviouralPredictions = simulated(behaviouralFilter,:);
        
    % Calculate the likelihoods
    display('   Calculatting the likelihoods and parameter distributions ...');
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
    for i=1: length(parameterNames)
        parameterDistributions.(parameterNames{i}) = sortrows([behaviouralParams(:,i),likelihood]);
        parameterDistributions.(parameterNames{i})(:,2) = cumsum(parameterDistributions.(parameterNames{i})(:,2));
    end
    
    
    display('   Calculating the simulation confidence intervals...');
    % calculate 90% confidence intervals for phosphorus predictions
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
                predictionCIs(predictionCIs(:,2)<effectiveZero,2) = -1e-9;
    end
    
    % Plot the simulated P
    figure()    
    plot(simDates, observedFilt(:,end),simDates, medianPrediction(:,2),simDates, predictionCIs(:,2),simDates, predictionCIs(:,3)); 
    legend('Observed','Modelled-Median', 'Modelled-5th', 'Modelled-95th');
    title(['GLUE Prediction uncertainty for ',num2str(realisations),' realsiation and a threshold of ', num2str(threshold)]);
    xlabel('Time (years)');
    if doMonthly
        ylabel('Streamflow (mm/month)');        
    else
        ylabel('Streamflow (mm/d)');        
    end
    box on;
    
    
    pcntObsBelowLowerCI = sum(observedFilt(:,end) < predictionCIs(:,2))/numberOutputTimes;
    pcntObsAboveUpperCI = sum(observedFilt(:,end) > predictionCIs(:,3))/numberOutputTimes;
    display('   Finished GLUE anaalysis.');
    display(['   The fraction of observations above the upper confidence interval =',num2str(pcntObsAboveUpperCI)]);
    display(['   The fraction of observations below the lower confidence interval =',num2str(pcntObsBelowLowerCI)]);
    
    simulated = simulated';
    behaviouralPredictions = behaviouralPredictions';
end

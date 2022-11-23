function [ performance, simulated, outDates ] = HBV_MonteCarlo( MCparameters, ...
    parameterStructure, startDate, endDate, climateData, observedFlow, doMonthly, keepTimeSeries, normaliseSSE)
%METHANEPRODUCTION_ Summary of this function goes here
%   Detailed explanation goes here

runsPerBlock = 10000;
parameterNames = fieldnames(parameterStructure);

% Irrespective of doMonthly, perform daily evaluations initially.  See
% below for later conversion to monthly if needed.
[~, ~, Q] = objectiveFunction( MCparameters(1,:), parameterNames,  ...
    parameterStructure, startDate, endDate, climateData, observedFlow, false);
numDays=length(Q(:,1));
dateCols = Q(:,1:3);

if doMonthly
    observedFlowMonthly = convertDailyToMonthly(observedFlow);
    Q = convertDailyToMonthly(Q);
end

outDates = Q(:,1:3);

% get dates
if doMonthly
    obsDates = datenum(observedFlowMonthly(:,1), observedFlowMonthly(:,2), observedFlowMonthly(:,3));
else
    obsDates = datenum(observedFlow(:,1), observedFlow(:,2), observedFlow(:,3));
end

modDates = datenum(outDates(:,1), outDates(:,2), outDates(:,3));

% Calculate the residuals
% Trim series to match start/end dates before calculating residuals
firstDate = max([modDates(1),obsDates(1), startDate]);
lastDate = min([modDates(end),obsDates(end), endDate]);
Qfilt = modDates >= firstDate & modDates <= lastDate;
obsFilt = obsDates >=firstDate & obsDates <= lastDate;

numTimes=length(Q(:,1));

numRealisations = length(MCparameters(:,1));
performance = NaN(numRealisations,1);

if keepTimeSeries
    simulated = NaN(numTimes,numRealisations);
else
    simulated = [];
end

% Run in block to save memory

blockBounds = runsPerBlock:runsPerBlock:numRealisations;
if isempty(blockBounds)
    blockBounds = numRealisations;
elseif blockBounds(end)<numRealisations
    blockBounds = [blockBounds, numRealisations];
end
blockBounds =[1 blockBounds(1:end-1)+1; blockBounds]';
numBlocks = length(blockBounds(:,1));

% Now run the blocks
for j=1:numBlocks
    lb = blockBounds(j,1);
    ub = blockBounds(j,2);
    
    simTemp = NaN(numDays,ub-lb+1);
    parTemp = MCparameters(lb:ub,:);

    parfor i=1:ub-lb+1
        parameterVector = parTemp(i,:);
        [sse, ~, Q] = objectiveFunction( parameterVector, parameterNames,  ...
            parameterStructure, startDate, endDate, climateData, observedFlow, false);
        simTemp(:,i)= Q(:,end);
        perfTemp(i) = sse;
    end
    
    if doMonthly
        % Need to recalculate sse and flow at monthly step as we did daily
        % evaluations above.  The reason for doing this is that the conversion
        % to  monthly is slow and can be done once only here, rather than for
        % every realisation.
        simTemp = convertDailyToMonthly([dateCols simTemp]);
        simTemp = simTemp(:,4:end);

        residuals = simTemp(Qfilt,:) - repmat(observedFlowMonthly(obsFilt,4),1,ub-lb+1);
        residuals = residuals(~isnan(residuals(:,1)),:);

        % Calculate the objective function
        perfTemp= residuals' * residuals;
        perfTemp = perfTemp(logical(eye(ub-lb+1)));

    end %doMonthly

	performance(lb:ub) = perfTemp;
    if keepTimeSeries
        simulated(:,lb:ub) = simTemp;
    end
end
if ~keepTimeSeries
    outDates = [];
end
% Normalise SSE if needed
if normaliseSSE
    if doMonthly
        observedDataFilter = ~isnan(observedFlowMonthly(:,end));
        meanObs = mean(observedFlowMonthly(observedDataFilter,end));
        observedSumSquares = (observedFlowMonthly(observedDataFilter,end)'-meanObs)* ...
            (observedFlowMonthly(observedDataFilter,end)-meanObs);
    else
        observedDataFilter = ~isnan(observedFlow(:,end));
        meanObs = mean(observedFlow(observedDataFilter,end));
        observedSumSquares = (observedFlow(observedDataFilter,end)'-meanObs)* ...
            (observedFlow(observedDataFilter,end)-meanObs);
    end
    performance = performance/observedSumSquares;
end

end


function [MCparameters, performance, simulated, outputDates, GLUE] = ...
    HBV_RegionalSensivityAnalysis(parameters, parametersLowerBound, ...
    parametersUpperBound, ClimateData, start_year, end_year, doMonthly, realisations, Observed_Flow, ...
    parameterNames, levels, varargin)
%
%HBV_REGIONALSENSITIVITYANALYSIS Regional Sensitivity Analysis of
%the HBV model.
%
%  This function undertaks a Regional Sensitivity Analysis (RSA) of the 
%  HBV model.
%
%  Runs are divided into ten groups based on the measures of performance 
%  supplied in performanceText.  Cumulative probability plots are produced
%  for each performance measure and parameter combination.  PlotMatrix 
%  style plots are also produced for the first and tenth decile groups for 
%  each performance measure.
%
%  [MCparameters, performance, simulated, outputDates, GLUE] =  ...
%       HBV_RegionalSensivityAnalysis(parameters, parametersLowerBound, ...
%       parametersUpperBound, ClimateData, start_year, end_year, ...
%       doMonthly, realisations, Observed_Flow, ...
%       parameterNames, levels)
%
%  INPUTS
%  parameters is a structure containing parameter values for the model.  It
%   is used to supply parameters for variables not included in the RSA
%
%  parametersLowerBound and parametersUpperBound define the limits between
%   which parameters are sampled from a uniform probability distribtuion.
%   Parameter sampling assumes independence between the parameters.
%   Parameter field names and field name order must be as found in
%   defaultParameterRanges.m
%
%  ClimateData is the climate data in the format required by HBV_noSnow.m
%
%  start_year, end_year specify the start and end times for the simulation
%
%  doMonthly controls whether a daily (false) or monthly (true) evaluation
%   of the objective function is used
%
%  realisations is the number of Monte Carlo runs undertaken for the RSA
%
%  Observed_Flow is a daily observed flow series in the format required by
%   HBV_noSnow.m
%
%  parameterNames is a cell matrix of parameter names corresponding to
%   fields in the various parameter structures that will be included in the
%   RSA
%
%  levels specifies the number of groups used in the RSA
%
%  OPTIONAL INPUTS ('name','value' pairs)
%   HBV_RegionalSensivityAnalysis(___,Name,Value)
%
% 'groupingType': 'linear' (default) or 'logarithmic'.  
%   Linear produces constant bin sizes, logarithmic has smaller bins for low
%   performance metric values - levels determinse the smallest bin size
%   Levels, Best runs bin size (% of runs)
%   2       64
%   3       32
%   4       16
%   5        8
%   6        4
%   7        2
%   8        1
%   9        0.5
%   10       0.25
%   etc
%
% 'GLUEbins': true, false (default)
%   If set true, GLUEthreshold is interpreted as the number of levels to
%   include in GLUE analysis
%
% 'GLUEthreshold': numerical
%   If a positive value is supplied, a GLUE analysis is conducted the and 
%   the threshold used in GLUE is the sum of squared errors divided by the
%   sum of squared variations of the observed from the the mean observed.  
%   i.e. 1 minus the Nash-Sutcliffe Coefficient of Efficiency
%   GLUE results are output are output in the structure GLUE
%
% 'keepTimeseries': true, false (default)
%   If true the modelled timeseries from the Regional are output.
%
% 'effectiveZero': numerical
%   If a positive value is supplied, any simulated values in a GLUE analysis 
%   are set to zero if they are less than effectiveZero
%
%
%  OUTPUTS
%  MCparameters is a matrix is parameters used in the model runs with each row
%   corresponding to a model run.  Dimensions are [realisations,7]
%
%  performance is a matrix of performance measures with dimensions
%   [realisations, length(performanceText]
%
%  simulated is a matrix of simulated flows for MCparameters
%
%  outputDates are the dates for the simulated flows
%
%  GLUE is a structure containting all the outputs from
%   HBV_GLUE_existingRuns.  These include:
%       medianPrediction
%       predictionCIs
%       parameterDistributions
%       behaviouralParams
%       behaviouralPredictions
%       likelihood
%       pcntObsAboveUpperCI
%       pcntObsBelowLowerCI
%   See GLUE_existingRuns for details
%
%  EXAMPLE
%
%  start_year = 1981;
%  end_year = 2006;
%  doMonthly = false;
%  realisations = 10000;
%  parameterNames = fieldnames(Params_initialModel);
%  levels = 10;
%  
%  [MCparameters, performance, simulated, outputDates, GLUE] =  ...
%     HBV_RegionalSensivityAnalysis(Params_initialModel, ...
%     Params_LowerBound, Params_UpperBound, sevensCreekClimate, start_year, end_year, ...
%     doMonthly, realisations, sevensCreek_mm_day, parameterNames, levels, ...
%     'groupingType','logarithmic', 'GLUEthreshold', 0.5, 'GLUEbins', false, ...
%     'keepTimeseries', false);
%
%
%  RELATED FUNCTIONS
%  parameterMatrices.m generates a matrix of parameters for input to HBV
%
%  HBV.m is the HBV model solver
%
%  HBV_GLUE_existingRuns undertakes GLUE analysis on a subset of the RSA
%  results
%
%  Written by Andrew Western, University of Melbourne, 29 April 2014

%% Start code

% Hardwire varaibles
% performanceText specifies the graph titles for each performance measture
% that is plotted
performanceText = {'norm SSE'};
performanceColumns = 1;

% work out boundaries for grouping and parse name value pairs
narginchk(11,21);

p = inputParser;
defaultGroupingType = 'linear';
validGrpTypes = {'linear','logarithmic'};
checkGrpType = @(x) any(validatestring(x,validGrpTypes));

defaultGLUEthreshold = -1;
defaultKeepTS = false;
defaultGLUEbins = false;
defaultEffectiveZero = 0;

addOptional(p,'groupingType',defaultGroupingType,checkGrpType);
addParameter(p,'GLUEthreshold',defaultGLUEthreshold,@isnumeric);
addParameter(p,'GLUEbins',defaultGLUEbins,@islogical);
addParameter(p,'keepTimeseries',defaultKeepTS,@islogical);
addParameter(p,'effectiveZero',defaultEffectiveZero,@isnumeric);


parse(p,varargin{:});

% round realisations to nearest factor of "levels" to get even group sizes
if strcmp(p.Results.groupingType,'linear')
    subsetLength=round(realisations/levels);
    realisations = subsetLength*levels;
    groupBounds = [(0:levels-1)*subsetLength+1;(1:levels)*subsetLength]';
    else
    groupBounds=2.^(1:levels-1);
    groupBounds = round(groupBounds(1)/groupBounds(end)*0.64*realisations)*2.^(0:levels-2); % results in bounds at 1%, 2% etc of runs
    groupBounds = [1, groupBounds+1;groupBounds, realisations]';
%     subsetLength=groupBounds(end,2)-groupBounds(end-1,2)+1;
    disp('HBV_RegionalSensivityAnalysis: Using logarithmic grouping');
    disp(['Run percentage bounds are:',sprintf(' %0.2g,',64*0.5.^(levels-2:-1:1)),' and 64% of runs']);
end

if p.Results.GLUEthreshold > 0
    doGLUE = true;
    threshold = p.Results.GLUEthreshold;
    if p.Results.GLUEbins
        GLUEusesRSAbins = true;
        threshold = round(threshold);
        disp(sprintf('GLUE will use the best %0.3g%% of runs',64*0.5^(levels-1-threshold)));
    else
        GLUEusesRSAbins = false;
    end
else
    doGLUE = false;
end

if p.Results.effectiveZero > 0
    forceToZero = true;
    effectiveZero = p.Results.effectiveZero;
else
    forceToZero = false;
    effectiveZero = 0;    
end
    
% extract parameter names from the parameter bound structure and find the
% columns in the parameter matrix that need to be replaced by stochastic
% parameters
allParamNames=fieldnames(parametersLowerBound);
numParams=length(allParamNames);
numRSAparams=length(parameterNames);
parameterColumns = NaN(1,numRSAparams);
for j=1:numRSAparams
    for i=1:numParams
        if strcmp(allParamNames(i),parameterNames(j)) 
            break
        end
    end
    parameterColumns(j)=i;
end

%%
%Setup the time series

% Trim to min of start climate and start obs flow data.
ClimateDates = datenum( ClimateData(:,1), ClimateData(:,2), ClimateData(:,3));
FlowDates = datenum( Observed_Flow(:,1), Observed_Flow(:,2), Observed_Flow(:,3));
startDate = max( max( min(ClimateDates) , min(FlowDates)), datenum(start_year, 1, 1) );
endDate = min( min( max(ClimateDates) , max(FlowDates)), datenum(end_year, 12, 31) );

FlowDates_filt = FlowDates >= startDate & FlowDates <= endDate;    
Observed_Flow_trimmed = Observed_Flow(FlowDates_filt,:);        

% set up parameter matrix
MCparameters = parameterMatrices(parameters, parametersLowerBound, ...
    parametersUpperBound, realisations, parameterNames);

% run the model and calculate performance measures for grouping runs
[performance, simulated, outputDates] = HBV_MonteCarlo(MCparameters, ...
    parameters, startDate, endDate, ClimateData, Observed_Flow, ...
    doMonthly, p.Results.keepTimeseries, true); 

% DRAW PLOTS
% First calculate an matrix of index values that code the numerical order
% of the performance measures.  Note that the order of the performance 
% matrix is not changed by this but the matrix indices records the
% numerical order from smallest to largest which is used to extract
% parameter sets from the params matrix
[~,indices]=sort(performance);

% Declare a matrix to store parameter values for CDF plotting purposes
% cumulative=NaN(subsetLength,levels,numRSAparams);

% Loop through the performance measures
for p=1:length(performanceText)
    % extract the parameter values associated with each decile of the
    % performance measure
    for i=levels:-1:1
%          cumulative(groupBounds(i,1):groupBounds(i,2),i,:)= ...
%              MCparameters(indices(groupBounds(i,1):groupBounds(i,2),performanceColumns(p)),parameterColumns);
        cumulative(i).data= sort(...
            MCparameters(indices(groupBounds(i,1):groupBounds(i,2),performanceColumns(p)),parameterColumns));
        probability(i).data=(1:1:length(cumulative(i).data(:,1)))/length(cumulative(i).data(:,1));
         
    end
    
    % sort the parameter values for plotting CDFs
%     cumulative=sort(cumulative);

    % extract the decile boundaries and set up labelling for the colour bar
    percentiles=[1;groupBounds(:,2)];
    percentiles=performance(indices(percentiles,p),p);
    for i=levels+1:-1:2
        colorBarLabels(i)={[sprintf('%3.1e %2g',percentiles(i),groupBounds(i-1,2)/realisations*100),'%']}; 
    end;
    colorBarLabels(1)={[sprintf('%3.1e %2g',percentiles(1),0),'%']};

    % calculate the cumlative probabilites to for CDF plotting purposes
%     probability=[1:1:subsetLength]/subsetLength;
    % Loope through each parameter and produce the CDF plot
    for j=1:numRSAparams
%         j=parameterColumns(i);
            theFigures(j)=figure();
            set(gca, 'ColorOrder', colormap(jet(levels)));
            hold all;
            for i=1:levels
%                 plot(cumulative(:,i,j),probability,...
%                     'LineWidth',2);
                plot(cumulative(i).data(:,j),probability(i).data,...
                    'LineWidth',2);
            end
            
            % add legend and various text labels
            h=colorbar('Ticks',(0:levels+1)/levels,'TickLabels',colorBarLabels);
            v=get(h,'title');
            set(v,'string',performanceText(p),'FontSize',14);
            xlabel(parameterNames(j),'FontSize',11);
            ylabel('Cumulative Proportion','FontSize',11);
            title('Regional Sensitivity Analysis','FontSize',14);
    end

    % draw plotmatrix style plot for upper decile
    theFigures(j+1)=figure();
    [Handle,AX,BigAx,P,PAx] = plotmatrix(MCparameters(indices(groupBounds(end,1):groupBounds(end,2),p),parameterColumns),'.');
    title(strcat('Highest group: ', performanceText(p)),'FontSize',14);
    % Loop through parameters and label the relevant subplot axes in the
    % plotmatrix
    for q=1:numRSAparams
        ylabel(AX(q,1),parameterNames(q),'FontSize',11);
        xlabel(AX(numRSAparams,q),parameterNames(q),'FontSize',11);
    end
    
    % draw plotmatrix style plot for lower decile
    theFigures(j+2)=figure();
    [Handle,AX,BigAx,P,PAx] = plotmatrix(MCparameters(indices(groupBounds(1,1):groupBounds(1,2),p),parameterColumns),'.');
    title(strcat('Lowest group: ', performanceText(p)),'FontSize',14);
    % Loop through parameters and label the relevant subplot axes in the
    % plotmatrix
    for q=1:numRSAparams
        ylabel(AX(q,1),parameterNames(q),'FontSize',11);
        xlabel(AX(numRSAparams,q),parameterNames(q),'FontSize',11);
    end
    
end % for p=1:length(performanceText)

%% Run GLUE if required

if doGLUE
    % filter the right parameters here.
    if GLUEusesRSAbins
        parameterValuesTemp = MCparameters(indices(groupBounds(1,1):groupBounds(threshold,2)),:);
        GLUEfilt = false(size(performance));
        GLUEfilt(indices(groupBounds(1,1):groupBounds(threshold,2)))=true;
    else
        GLUEfilt = performance <= threshold;
        parameterValuesTemp = MCparameters(GLUEfilt,:);
    end %if GLUEusesRSAbins
    
    if sum(GLUEfilt)==0
        warning('There were no behavioural GLUE simulations - trying increasing the threshold');
        GLUE=[];
    else
        [tempPerformance, tempSimulated, theDates] = HBV_MonteCarlo(parameterValuesTemp, ...
            parameters, startDate, endDate, ClimateData, Observed_Flow_trimmed, doMonthly, true, true); 

        if doMonthly
            Observed_Flow_trimmed = convertDailyToMonthly(Observed_Flow_trimmed);
        end

        if GLUEusesRSAbins
            GLUE.threshold = max(tempPerformance);
        else
            GLUE.threshold = threshold;
        end

        [ GLUE.medianPrediction, GLUE.predictionCIs, GLUE.parameterDistributions, GLUE.behaviouralParams, ...
            GLUE.behaviouralPredictions, GLUE.likelihood, GLUE.pcntObsAboveUpperCI, GLUE.pcntObsBelowLowerCI ] = ...
            GLUE_existingRuns( parameterNames, parameterValuesTemp, ...
            tempPerformance, tempSimulated, GLUE.threshold, Observed_Flow_trimmed, ...
            doMonthly, forceToZero, effectiveZero);

        GLUE.outputDates = theDates;
    end %if sum(GLUEfilt)==0
    
    
else
    GLUE = [];
end %if doGLUE

end


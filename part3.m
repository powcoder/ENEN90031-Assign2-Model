% clear all
% close all
% load('assignment2_Data.mat')

%% Part 3

realisations = [];
doMonthly = false;

parameterNames =  {'fc';'beta';'pwp';'l'; 'k0'; 'k1'; 'kp';'k2'};

timeStart = datenum(2005,1,1);
timeEnd = datenum(2018,12,31);
start_year = year(timeStart);
end_year = year(timeEnd);

% Add from here to run GLUE

%Save key outputs from GLUE - you need to use the actual variable names you
%have in your code!
save('glue_output.mat','MCparameters', 'performance', 'simulated', ...
    'outDates', 'parameterNames', 'medianPrediction', 'predictionCIs', ...
    'parameterDistributions', 'behaviouralParams', ...
    'behaviouralPredictions', 'likelihood', 'pcntObsAboveUpperCI', 'pcntObsBelowLowerCI');


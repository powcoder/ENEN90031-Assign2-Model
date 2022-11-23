% clear all
% close all
% load('assignment2_Data.mat')

%% Part 2

realisations = [];
doMonthly = false;

timeStart = datenum(2005,1,1);
timeEnd = datenum(2018,12,31);

start_year = year(timeStart);
end_year = year(timeEnd);

parameterNames =  {'fc';'beta';'pwp';'l'; 'k0'; 'k1'; 'kp';'k2'};

% Add from here to run regional sensitivity analysis




% Save your Monte Carlo run outputs - you may need to alter variable names
% to those you used.
save('monteCarloRuns.mat', 'MCparameters', 'performance', 'simulated', 'outDates');

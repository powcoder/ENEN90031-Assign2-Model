%% Part 4
% clear all
% close all
% load('assignment2_Data.mat')

doMonthly = false;

timeStart = datenum(1988,1,1);
timeEnd = datenum(2018,12,31);
start_year = year(timeStart);
end_year = year(timeEnd);

head = 300;
catchmentArea = 361;

budget= []; % I suggest you start with broad budget range say 1e6 to 1e9 and modify

warning('off','all'); % this suppresses warnings from HBV_MonteCarlo regarding spinup

% Run hydroDesign

% Create climate change inputs

% Run hydroPower for climate change case


warning('on','all');
clear all
close all
load('assignment2_Data.mat')

doMonthly = false;
parameterNames =  {'fc';'beta';'pwp';'l'; 'k0'; 'k1'; 'kp';'k2'};

% Set Optimisation Settings
optimisationSettings.maxn = 200000;
optimisationSettings.kstop = 5;    
optimisationSettings.pcento = 1e-6;    
optimisationSettings.peps = 1e-6;
optimisationSettings.ngs = []; % You need to set this value

% Set dates
timeStart = datenum(2005,1,1);
timeEnd = datenum(2018,12,31);

% Add from here to run SCE
function [Storage, reliabitityTime, reliabilityDemand, Power] = hydroPower(streamFlowDates, ...
    streamFlow, demand, head, catchmentArea, storageCapacity, hydroCapacity, S0)

% HYDROPOWER is a function used to estimate the hydro power generated for
% an assumed capacity and to find the reliability of the dam
%
% SYNTAX
% [Storage, reliabitityTime, reliabilityDemand] = hydroPower(streamFlow,head,catchmentArea...
%   ,demand, capacity)
%
% INPUTS
% streamFlowDates: rows with year, month, day
% streamFlow : daily streamflow in (mm/day).  Multiple columns represent
%   multiple realisations
% demand: rows with year, month, day, demand (Watts)
% % head : elevation difference between dam and hydro power plant, (m)
% catchmentArea : total area of the catchment contributing to the
%                 streamflow , (km^2)
% capacity : dam capacity in volume (ML)
% hydroCapacity : maximum rated hydroelecticity output (W)
% S0 : initial storage in dam at the start of simulation (ML)
%
% OUTPUT
% Storage : daily storage in dam (ML)
% reliabityTime : calculated as, sum (days when powerGenerated >
% demand)/total no of days (-)
% reliabilityDemand : calculated as, 
%   sum (max(demand-powerGenerated),0)/sum(demand) (-)
%
% ASSUMPTIONS
% Evaporative loss from the dam is assumed to be zero.
% Change in head of water stored in the dam is assumed to be negligible
% compared to the elevation difference between dam and hydro electric
% station. i.e., this model is written for constant head condition.
% No releases are made other than for power generation or spills when the
% dam is full

% Check dates match
flowDates = datenum(streamFlowDates(:,1),streamFlowDates(:,2),streamFlowDates(:,3));
demandDates = datenum(demand(:,1),demand(:,2),demand(:,3));

startDate = max(flowDates(1),demandDates(1));
endDate = min(flowDates(end),demandDates(end));
flowFilt = flowDates >= startDate & flowDates<=endDate;
demandFilt = demandDates >= startDate & demandDates<= endDate;

flowDates = flowDates(flowFilt);
streamFlow = streamFlow(flowFilt,:);
demandDates = demandDates(demandFilt);
demand = demand(demandFilt,:);


try 
    dateDifference = flowDates - demandDates;
catch
    error('hydroPower: problem with flow and demand date matching');
end

if any(abs(dateDifference)>0.001);
    error('hydroPower: problem with flow and demand date matching');
end

% initialization
Storage = nan(size(streamFlow));
Vin=zeros(size(streamFlow)); %  [m^3] volume that flows into the dam
Vout=zeros(size(streamFlow)); %  [m^3] volume that flows out the dam


gamma = 9810; % [kg/m^2/s^2] density of water x gravity
% Q_turbine=20; % [m^3/s] flow per 1 hydrounit <<< should be changed according to the catchment char>>
E =0.92; % assuming 92% efficiency
% Power equation
% Power = E*gamma*Q*head; %(Watt)
timestep=3600*24; % one day of seconds
VoutMax = hydroCapacity/(E*gamma*head)*timestep / 1000;  %ML


Power = nan(size(streamFlow));
temp1 = nan(size(streamFlow));
temp = zeros(1,size(streamFlow,2));

numSteps = size(streamFlow,1);

for i=1:numSteps
    Vin(i,:)=streamFlow(i,:)*catchmentArea; % inflow volume (ML)
    VolOut=(demand(i,4)/(gamma*head*E))*timestep/1000; %%  Power equation : Power = E*gamma*Q*head; %(Watt)
    Vout(i,:)=min(min(VolOut,(S0+Vin(i,:))), VoutMax); % outflow volume limited 
        %to available water and available power plant capacity
        
    Storage(i,:) = min(storageCapacity,(S0+Vin(i,:)-Vout(i,:))); % dam storage at i th timestep
    S0=Storage(i,:);
    
    % reliability calculation
    Power(i,:) = E*gamma*head*(Vout(i,:)*1000)/timestep;
    powerFilt = Power(i,:)+1e-6<=demand(i,4); % these are failures to meet demand, the 1e-6 is to deal with rounding error.
    temp(powerFilt)=temp(powerFilt)+1;
    temp1(i,powerFilt)=demand(i,4)-Power(i,powerFilt);
    temp1(i,~powerFilt)=0;

end

reliabitityTime = 1-temp/numSteps;
reliabilityDemand =1-sum(temp1)./sum(demand(:,4));

end

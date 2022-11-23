function [bestReliability,bestHydroCapacity,bestDamCapacity] = ...
    hydroDesign(budget, streamFlowDates, streamFlow, demand, head, catchmentArea)
%HYDRODESIGN determines the division of capital between storage and
% other infrastructure that maximises reliability subject to a budget
% constraint.
%   Given a series of inflows and electrical demands, hydro electric
%   generation and storage of water in a reservoir are simulated to
%   determine reliability.  Based on a budge, investment is divided between
%   water storage infrastructure and other infrastructure in such a way
%   that reliability, defined as energy supplied divdided by energy
%   demanded, is maximised.  It searches for the best design using a grid
%   search.  The best design is the one that maximises reliability averaged
%   across all streamflow realisations.
%
%   SYNTAX
%   [bestReliability,bestHydroCapacity,bestDamCapacity] = ...
%       hydroDesign(budget, streamFlowDates, streamFlow, demand, 
%       head, catchmentArea)
%
%   INPUTS
%
%   budget: available budget ($)
%   streamFlowDates: a three column matrix with year, month, day columns
%       snd n rows
%   streamFlow: an n by m matrix where columns represent realisations of
%   dam inflow (mm)
%   demand: an n*4 matrix with year, month, day, electricity demand (W)
%   head: generating head across turbine (m)
%   catchmentArea: catchment area for inflow (km2)
%
%   OUTPUTS
%
%   bestReliability: a 1 by m vector of system reliability for the selected
%       design
%   bestHydroCapacity: the selected maximum generating capacity of the 
%       power station (W)
%   bestDamCapacity: the selected reservoir capacity (ML)
%
%   Written by Andrew Western, Univeristy of Melbourne, May 2020

% Set up controlling variables.
numCases = length(streamFlow(1,:));
searchStep = 0.1;
pattern = (-5:5)';
bestFraction = 0.5;
numRows=length(pattern);

% Declare search output
reliabilityDemand = nan(numRows,numCases);


%% search loop
for jj=1:3
    
    % grid of budget fraction
    % Force it to [0, 1] to avoid truncation errors taking it outside valid
    % range
    damBudgetFraction = min(1,max(0,searchStep*pattern + bestFraction)); 

    % divide budget and calculate size of dam and power staation
    
    % storage cost  $AU/ML = 13138vol^-0.555.   Vol in GL  
    % Petheram and McMahon https://doi.org/10.1016/j.hydroa.2019.100026
    % or $AU/ML = 607476 Vol ^ -0.555.   Vol in ML
    damCapacity = (budget*damBudgetFraction / 607476).^(1/(1-0.555));  %ML

    % 3000  $/KW  (% Large US is around $2500 / KW excl. Dam for a 500MW development.
    % IRENA, Cost of Hydro 2012
    hydroCapacity = budget*(1-damBudgetFraction) / 3; %W

    % Calculate grid of reliability for various budget allocations
    for ii=1:numRows
        [,~,reliabilityDemand(ii,:),~] = hydroPower(streamFlowDates, ...
            streamFlow, demand, head, catchmentArea, damCapacity(ii,:), ...
            hydroCapacity(ii,:), damCapacity(ii,:)/2);
    end
    
    % Find the budget allocation with best performance and refine search
    % grid to include intervals either side of the best.
    
    %Reliability average across realisations
    meanReliabilityDemand = mean(reliabilityDemand,2); 
    % find maximum reliability realisations 
    maxFilt = meanReliabilityDemand == max(meanReliabilityDemand); 
    r=find(maxFilt,1);
    numMax = sum(maxFilt); % check how many maximimum values there are
    if numMax > 1 
        %more than one maxima, find maxima closest the the average location
        %of all maxima
        r=floor(mean(find(maxFilt(:))));
        if ~maxFilt(r)
            for ll=1:5
                if any(maxFilt(max(r-ll,1):min(r+ll,numRows)))
                    r = r-ll+find(maxFilt(max(r-ll,1):min(r+ll,numRows)),1)-1;
                    break
                end
            end
        end
        if jj<3
            % make sure it is not the first or last point, except for last
            % loop
            r = max(2,min(10,r)); 
        end
    else
        if jj<3
            % make sure it is not the first or last point, except for last
            % loop - the search pattern extends from point before to point
            % after r
            r=max(2,min(10,find(maxFilt)));
        end
    end
    % Reset the maxFilt to reflect position determined
    maxFilt(:) = false;
    maxFilt(r) = true;
    
    % Update best fraction etc
    bestFraction = damBudgetFraction(maxFilt,:);
    bestDamCapacity = damCapacity(maxFilt);
    bestHydroCapacity = hydroCapacity(maxFilt);
    bestReliability = reliabilityDemand(maxFilt,:);
    searchStep=searchStep/5;
end


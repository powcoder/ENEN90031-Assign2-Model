%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bestx,bestf] = sceua(funcHandle, x0,bl,bu,maxn,kstop,pcento,peps,ngs,iseed,iniflg, varargin)

% This is the subroutine implementing the SCE algorithm, 
% written by Q.Duan, 9/2004
%
% Modified by Tim Peterson to include the following features:
%	- input of a function handle, rather than always using the objective function called 'functn.m'.
%	- input of variable arguments for passing data to the objective function
%	- parralelisation of the simplex step of each complex.
%
% Definition:
%  funcHandle = function handle to the model to be ran.
%  x0 = the initial parameter array at the start;
%     = the optimized parameter array at the end;
%  f0 = the objective function value corresponding to the initial parameters
%     = the objective function value corresponding to the optimized parameters
%  bl = the lower bound of the parameters
%  bu = the upper bound of the parameters
%  iseed = the random seed number (for repetetive testing purpose)
%  iniflg = flag for initial parameter array (=1, included it in initial
%           population; otherwise, not included)
%  ngs = number of complexes (sub-populations)
%  npg = number of members in a complex 
%  nps = number of members in a simplex
%  nspl = number of evolution steps for each complex before shuffling
%  mings = minimum number of complexes required during the optimization process
%  maxn = maximum number of function evaluations allowed during optimization
%  kstop = maximum number of evolution loops before convergency
%  percento = the percentage change allowed in kstop loops before convergency
%  varargin = cell strcuture to allow input of model specific data, eg forcing data or obs. data 

% LIST OF LOCAL VARIABLES
%    x(.,.) = coordinates of points in the population
%    xf(.) = function values of x(.,.)
%    xx(.) = coordinates of a single point in x
%    cx(.,.) = coordinates of points in a complex
%    cf(.) = function values of cx(.,.)
%    s(.,.) = coordinates of points in the current simplex
%    sf(.) = function values of s(.,.)
%    bestx(.) = best point at current shuffling loop
%    bestf = function value of bestx(.)
%    worstx(.) = worst point at current shuffling loop
%    worstf = function value of worstx(.)
%    xnstd(.) = standard deviation of parameters in the population
%    gnrng = normalized geometri%mean of parameter ranges
%    lcs(.) = indices locating position of s(.,.) in x(.,.)
%    bound(.) = bound on ith variable being optimized
%    ngs1 = number of complexes in current population
%    ngs2 = number of complexes in last population
%    iseed1 = current random seed
%    criter(.) = vector containing the best criterion values of the last
%                10 shuffling loops


% Output detailed messages to the user.
beVerbose = true;

% Initialize SCE parameters:
nopt=size(x0,2);
npg=2*nopt+1;
nps=nopt+1;
nspl=npg;
mings=ngs;
npt=npg*ngs;

bound = bu-bl;

% Create an initial population to fill array x(npt,nopt):
% If x0 contains npt parameter sets, then don't re-sample initial parameter
% sets
disp(['Number of points per simplex = ', num2str(nps)]);
disp(['Number of complexes = ', num2str(ngs)]);
disp(['Total number of random parameter sets = ', num2str(npt)]);
disp('Randomly sampling the initial parameters...');
rand('seed',iseed);
if iniflg==1 && size(x0,1) ==1 
    x(1,:)=x0; 
elseif iniflg==1 && size(x0,1) >=npt 
    x = x0(1:npt,:);
else
    x=zeros(npt,nopt);
    parfor i=1:npt;
        x(i,:)=bl+rand(1,nopt).*bound;
    end;
end


nloop=0;
icall=0;
xf = inf(1,npt);
icall_temp = zeros(1,npt);
disp('Calculating the objective function value at each random parameter set...');
parfor i=1:npt;
    nModelFails = 0;
    while 1
        try
            icall_temp(i) = icall_temp(i) + 1;
            xf(i) = feval(funcHandle,x(i,:), varargin{:});
            
            % Resample parameters if Inf or NaN
            if isinf(xf(i)) || isnan(xf(i))
                x(i,:)=bl+rand(1,nopt).*bound;
            else
                break;
            end
        catch modelError
            % Resample parameter set, then re-run model.
            x(i,:)=bl+rand(1,nopt).*bound;
            nModelFails = nModelFails + 1;
        end
        
        if nModelFails>100;
            error('The model failed for 100 parameter sets. Check the objective function, model and parameters.');
        end
    end
    
end;
f0=xf(1);
disp('Finsihed calculating the objective function value at each random parameter set.');
disp('Starting SCE calibration...');

% Sum function calls
icall = sum(icall_temp);
clear icall_temp

% Sort the population in order of increasing function values;
[xf,idx]=sort(xf);
x=x(idx,:);

% Record the best and worst points;
bestx=x(1,:); bestf=xf(1);
worstx=x(npt,:); worstf=xf(npt);

% Compute the standard deviation for each parameter
xnstd=std(x);

% Computes the normalized geometric range of the parameters
gnrng=exp(mean(log((max(x)-min(x))./bound)));

if beVerbose
	disp('The SCE Initial Loop: 0');
	disp(['BESTF  : ' num2str(bestf)]);
	disp(['BESTX  : ' num2str(bestx)]);
	disp(['WORSTF : ' num2str(worstf)]);
	disp(['WORSTX : ' num2str(worstx)]);
	disp(' ');

	% Check for convergency;
	if icall >= maxn;
	    disp('*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE LIMIT');
	    disp('ON THE MAXIMUM NUMBER OF TRIALS ');
	    disp(maxn);
	    disp('HAS BEEN EXCEEDED.  SEARCH WAS STOPPED AT TRIAL NUMBER:');
	    disp(icall);
	    disp('OF THE INITIAL LOOP!');
	end;

	if gnrng < peps;
	    disp('THE POPULATION HAS CONVERGED TO A PRESPECIFIED SMALL PARAMETER SPACE');
	end;
else
    % Output summary of intial results
    disp(['        - SCE Evolution Loop ', num2str(0), ' Objective Function Results: - Best: ', num2str(bestf), ' - Worst: ', num2str(worstf), ' - Inter-quantile range: ', num2str(iqr(xf))]);   
end
% Begin evolution loops:
nloop = 0;
criter=[];
criter_change=1e+5;

while icall<maxn & gnrng>peps & criter_change>pcento;
    nloop=nloop+1;
    
    % Initialise function calls per complex
    complex_fcalls = zeros(ngs,1);

    % Partition the population into complexes (sub-populations);
    cx = zeros(npg, nopt, ngs);
    cf = zeros(nopt, ngs);
    for igs = 1: ngs        
        k1=1:npg;
        k2=(k1-1)*ngs+igs;
        cx(k1,:,igs) = x(k2,:);
        cf(k1,igs) = xf(k2);        
    end
    
    % Loop on complexes (sub-populations);
    parfor igs = 1: ngs;
        [cx_temp{igs}, cf_temp{igs},  complex_fcalls(igs)] = evolve_complex(nspl, npg, nps, cx(:,:,igs), cf(:,igs) , funcHandle, bl,bu,maxn, varargin{:});
    end;
    
    % Replace the complex back into the population;    
    for igs = 1: ngs;        
        k1=1:npg;
        k2=(k1-1)*ngs+igs;              
        
        x(k2,:) = cx_temp{igs};
        xf(k2) = cf_temp{igs};                
    end
    
    % Sum of function calls and add to icall
    icall = icall + sum(complex_fcalls);
    
    % Shuffled the complexes;
    [xf,idx] = sort(xf); x=x(idx,:);
    %PX=x; PF=xf;
    
    % Record the best and worst points;
    bestx=x(1,:); bestf=xf(1);
    worstx=x(npt,:); worstf=xf(npt);
    %BESTX=[BESTX;bestx]; BESTF=[BESTF;bestf];ICALL=[ICALL;icall];

    % Compute the standard deviation for each parameter
    xnstd=std(x);

    % Computes the normalized geometric range of the parameters
    gnrng=exp(mean(log((max(x)-min(x))./bound)));

    if beVerbose
        disp(['Evolution Loop: ' num2str(nloop) '  - Trial - ' num2str(icall)]);
        disp(['BESTF  : ' num2str(bestf)]);
        disp(['BESTX  : ' num2str(bestx)]);
        disp(['WORSTF : ' num2str(worstf)]);
        disp(['WORSTX : ' num2str(worstx)]);
        disp(['Interquartile range of F : ' num2str(iqr(xf))]);
        disp(' ');
    else
        disp(['        - SCE Evolution Loop ', num2str(nloop), ' Objective Function Results: - Best: ', num2str(bestf), ' - Worst: ', num2str(worstf), ' - Inter-quantile range: ', num2str(iqr(xf))]);
    end

    % Check for convergency;
    if beVerbose
        if icall >= maxn;
            disp('*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE LIMIT');
            disp(['ON THE MAXIMUM NUMBER OF TRIALS ' num2str(maxn) ' HAS BEEN EXCEEDED!']);
        end;

        if gnrng < peps;
            disp('THE POPULATION HAS CONVERGED TO A PRESPECIFIED SMALL PARAMETER SPACE');
        end;
    end

    criter=[criter;bestf];
    if (nloop >= kstop);
        criter_change=abs(criter(nloop)-criter(nloop-kstop+1))*100;
        criter_change=criter_change/mean(abs(criter(nloop-kstop+1:nloop)));
        if criter_change < pcento & beVerbose;
            disp(['THE BEST POINT HAS IMPROVED IN LAST ' num2str(kstop) ' LOOPS BY ', ...
                  'LESS THAN THE THRESHOLD ' num2str(pcento) '%']);
            disp('CONVERGENCY HAS ACHIEVED BASED ON OBJECTIVE FUNCTION CRITERIA!!!')
        end;
    end;
    
% End of the Outer Loops
end;

if beVerbose
    disp(['SEARCH WAS STOPPED AT TRIAL NUMBER: ' num2str(icall)]);
    disp(['NORMALIZED GEOMETRIC RANGE = ' num2str(gnrng)]);
    disp(['THE BEST POINT HAS IMPROVED IN LAST ' num2str(kstop) ' LOOPS BY ', ...
           num2str(criter_change) '%']);
end
% END of Subroutine sceua
return;
end

function [cx, cf, complex_fcalls] = evolve_complex(nspl, npg, nps, cx, cf, funcHandle, bl,bu,maxn, varargin)

        % Evolve sub-population igs for nspl steps:
        complex_fcalls = 0;
        for loop=1:nspl;
            
            % Select simplex by sampling the complex according to a linear
            % probability distribution
            lcs(1) = 1;
            for k3=2:nps;
                for iter=1:1000;
                    lpos = 1 + floor(npg+0.5-sqrt((npg+0.5)^2 - npg*(npg+1)*rand));
                    idx=find(lcs(1:k3-1)==lpos); 
                    if isempty(idx); 
                        break; 
                    end
                end
                lcs(k3) = lpos;
            end
            lcs=sort(lcs);

            % Construct the simplex:
            %s = zeros(nps,nopt);
            s=cx(lcs,:); 
            sf = cf(lcs);    
            
            [snew,fnew, fcalls]=cceua(funcHandle, s,sf,bl,bu,maxn, varargin{:});

            % Increase function call count.
            complex_fcalls = complex_fcalls + fcalls;
            
            % Replace the worst point in Simplex with the new point:
            s(nps,:) = snew; 
            sf(nps) = fnew;            
            
            % Replace the simplex into the complex;
            cx(lcs,:) = s;
            cf(lcs) = sf;
            
            % Sort the complex;
            [cf,idx] = sort(cf); 
            cx=cx(idx,:);
            
        % End of Inner Loop for Competitive Evolution of Simplexes
        end;

end

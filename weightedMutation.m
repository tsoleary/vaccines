function mutationChildren = weightedMutation(parents,options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation,mutationRate, weight)
%MUTATIONUNIFORM Uniform multi-point mutation.
%   MUTATIONCHILDREN = MUTATIONUNIFORM(PARENTS,OPTIONS,GENOMELENGTH,...
%                      FITNESSFCN,STATE,THISSCORE,THISPOPULATION, ...
%                      MUTATIONRATE) Creates the mutated children using
%   uniform mutations at multiple points. Mutated genes are uniformly 
%   distributed over the range of the gene. The new value is NOT a function
%   of the parents value for the gene.
%
%   Example:
%     options = optimoptions('ga','MutationFcn', @mutationuniform); 
%
%   This will create an options structure specifying the mutation
%   function to be used is MUTATIONUNIFORM.  Since the MUTATIONRATE is
%   not specified, the default value of 0.01 is used.
%
%     mutationRate = 0.05;
%     options = optimoptions('ga','MutationFcn', {@mutationuniform, mutationRate});
%
%   This will create an options structure specifying the mutation
%   function to be used is MUTATIONUNIFORM and the MUTATIONRATE is
%   user specified to be 0.05.
%
%   (Note: If calling gamultiobj, replace 'ga' with 'gamultiobj' in the
%   above examples)

%   Copyright 2003-2015 The MathWorks, Inc.
%   % EDITED ON NOV 10, 2019 BY ALEX BURNHAM TO FAVOR 0 MUTATIONS WITH A
%     GIVEN PROBABILTY: WEIGHT 

if nargin < 8 || isempty(mutationRate)
    mutationRate = 0.01; % default mutation rate
end

if(strcmpi(options.PopulationType,'doubleVector'))
    mutationChildren = zeros(length(parents),GenomeLength);
    for i=1:length(parents)
        child = thisPopulation(parents(i),:);
        % Each element of the genome has mutationRate chance of being mutated.
        mutationPoints = find(rand(1,length(child)) < mutationRate);
        % each gene is replaced with a value chosen randomly from the range.
        range = options.PopInitRange;
        % range can have one column or one for each gene.
        [r,c] = size(range);
        if(c ~= 1)
            range = range(:,mutationPoints);
        end   
        lower = range(1,:);
        upper = range(2,:);
        span = upper - lower;
        
        
        child(mutationPoints) = lower + rand(1,length(mutationPoints)) .* span;
        mutationChildren(i,:) = child;
    end
elseif(strcmpi(options.PopulationType,'bitString'))
    
    mutationChildren = zeros(length(parents),GenomeLength);
    for i=1:length(parents)
        child = thisPopulation(parents(i),:);
        mutationPoints = find(rand(1,length(child)) < mutationRate);
        child(mutationPoints) = randsample([0 1], 1, true, [weight 1-weight]);  % this is where mutation happens (opposite of what is there) (previous: ~child(mutationPoints);)
        mutationChildren(i,:) = child;
    end
    
end

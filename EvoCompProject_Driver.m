% Evolutionary Computation
% Project Driver
% 10, November 2019
% Blake Williams, Thomas O'Leary and Alex Burnham

% Starting parameters:
P = 5;
N = 10;
mutProb = 1/N; % probabilty of mutation
probZero = 0.7; % probability allele is 0
thresh = 3; % set our threashold for what our weighted genome sums to

% create weighted binary genome
genomeFlat = randsample([0 1], P*N, true, [probZero 1-probZero]);
genomeMat = vec2mat(genomeFlat, N); % creates weighted genome matrix with N cols

% run the normalizeation function (this goes in Blake's fitness function)
normGenome = NormalizeGenome(genomeMat, thresh); 




% call fitness function and adj_list with anonoymous function 
% NETWORK is still required (something from Blake?)
vaccineFitness = @(pop)(FitnessFunction(pop, adj_list_network));

% set options for ga toolbox for bitstring
vaccineOpts = gaoptimset(...
    'PopulationType', 'bitstring', ...
    'InitialPopulation', genomeMat, ...
    'CrossoverFcn', @crossoversinglepoint, ...
    'CrossoverFraction', 0.2, ...
    'SelectionFcn',{@selectiontournament,2}, ...
    'Generations', 500, ...
    'Vectorized', 'on', ...
    'FitnessLimit', 0, ...
    'MutationFcn', {@weightedMutation mutProb});

% run GA (fitnessFunction, genomeLength, passOptions)
[x, fval] = ga(vaccineFitness, N, vaccineOpts);










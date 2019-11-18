% Evolutionary Computation
% Project Driver
% 10, November 2019
% Blake Williams, Thomas O'Leary and Alex Burnham

% Starting parameters:

adj_mat_network=A; %%

P = 10;
N = 100;
mutProb = 1/N; % probabilty of mutation
probZero = 0.9; % probability allele is 0
thresh = 3; % set our threashold for what our weighted genome sums to
threshold = .5;
transcendence = 2;

% create weighted binary genome
genomeFlat = randsample([0 1], P*N, true, [probZero 1-probZero]);
genomeMat = vec2mat(genomeFlat, N); % creates weighted genome matrix with N cols

% run the normalizeation function (this goes in Blake's fitness function)
%normGenome = NormalizeGenome(genomeMat, thresh); 




% call fitness function and adj_list with anonoymous function 
% NETWORK is still required (something from Blake?)
%vaccineFitness = @(pop)(SpreadingFitnessFcn(pop, adj_mat_network, threshold, transcendence ));
 


% set options for ga toolbox for bitstring
vaccineOpts = gaoptimset(...
    'PopulationType', 'bitstring', ...
    'InitialPopulation', genomeMat, ...
    'CrossoverFcn', @crossoversinglepoint, ...
    'CrossoverFraction', 0.2, ...
    'SelectionFcn',{@selectiontournament,2}, ...
    'Vectorized','on',...
    'Generations', 500, ...
    'FitnessLimit', 0, ...
    'MutationFcn', {@mutationuniform, mutProb});
%s{@weightedMutation...
                   %  'mutatinRate',mutProb,...
                   %  'weight', probZero}
    

% run GA (fitnessFunction, genomeLength, passOptions)
%[x, fval] = ga(vaccineFitness, N, vaccineOpts);

[x,fval] = ga(@SpreadingFitnessFcn, N, vaccineOpts);

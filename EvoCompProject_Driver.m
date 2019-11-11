% Evolutionary Computation
% Project Driver
% 10, November 2019
% Blake Williams, Thomas O'Leary and Alex Burnham

% Starting parameters:
P = 5;
N = 10;
m = 10;
mutProb = 1/N; % probabilty of mutation
parents = randi(P, [1,m]); % indices m less than or equal to N and unique rows
probZero = 0.7; % probability allele is 0
thresh = 3; % set our threashold for what our weighted genome sums to

% create weighted binary genome
genomeFlat = randsample([0 1], P*N, true, [probZero 1-probZero]);
genomeMat = vec2mat(genomeFlat, N); % creates weighted genome matrix with N cols

% run the normalizeation function 
normGenome = NormalizeGenome(genomeMat, thresh); 

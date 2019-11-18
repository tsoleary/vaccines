function [supercritical_P]=SpreadingFitnessFcn(adjacency_mat, vaccine_vector, threshold, transcendence)
    % PARAMETERS
    % adjacency_mat      : adjmat of the genotype network
    % vaccine_vector     : length N= # strains, in [0 1]. 0 if nobody is vaccinated, 1 if fully vaccinated
    % epidemic threshold : from 0 to 1, should be ~0.5 (or 1/[transmission rate of
    %                      some disease]).

    % # nodes
    N=size(adjacency_mat,1);

    %strain transmission coefficients
    W=ones(1,N);

    %distance between all nodes (shortest paths)
    D=distances(graph(adjacency_mat));

    %compute reduction in betas due to transcendence across all nodes
    transcending_immunity= 1- exp(-D(find(vaccine_vector),:)/transcendence)';
    
    %weight effects by vaccine strength (from 0 to 1 for each vaccine
    % strain).
    transcending_immunity= 1- vaccine_vector(find(vaccine_vector)) .* (1-transcending_immunity);
    
    %take product for overlapping effects of multiple vaccination strains
    W=prod(transcending_immunity,2)';

    % THE FITNESS:
    %find proportion of supercritical strains: minimize this proportion
    supercritical_P=(sum(W>threshold))/N;
    
    

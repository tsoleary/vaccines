function [supercritical_P]=SpreadingFitnessFcn(adjacency_mat, vaccine_vector, vaccine_p, threshold, transcendence)
    % PARAMETERS
    % adjacency_mat      : adjmat of the genotype network
    % vaccine_vector     : length N, 0 if not vaccinated, 1 if vaccinated
    % vaccine_p          : from 0 to 1,   1:= transmission->0, 
    %                                   0.5:= transmission->1/2.
    % (epidemic) threshold : from 0 to 1, should be 0.5 (or 1/[transmission rate of
    %                      some disease]).
    % transcendence : from 0 to 10+, try 1 to 10 to start. Higher= wider spread of immunity, lower=more localized.


    % # nodes
    N=size(adjacency_mat,1);

    %strain transmission coefficients
    W=ones(1,N);

    %distance between all nodes (shortest paths)
    D=distances(graph(adjacency_mat));

    %compute reduction in betas due to transcendence across all nodes
    transcending_immunity= 1- vaccine_p * exp(-D(find(vaccine_vector),:)/transcendence)';
    
    %take product for overlapping effects of multiple vaccination strains
    W=prod(transcending_immunity,2)';

    % THE FITNESS:
    %find proportion of supercritical strains: minimize this proportion
    supercritical_P=(sum(W>threshold))/N;
    

function mean_compsize_P=SpreadingFitnessFcnCompSizeGrow(vaccine_pop, adjacency_mat, full_adjacency_mat, threshold, transcendence)
    % PARAMETERS
    % adjacency_mat      : adjmat of the genotype network
    % vaccine_vector     : length N= # strains, in [0 1]. 0 if nobody is vaccinated, 1 if fully vaccinated
    % epidemic threshold : from 0 to 1, should be ~0.5 (or 1/[transmission rate of
    %                      some disease]).
    mean_compsize_P=zeros(1,size(vaccine_pop,1));
    for loop=1:size(vaccine_pop,1)
        vaccine_vector=vaccine_pop(loop,:);
        % # nodes
        N=size(adjacency_mat,1);

        %strain transmission coefficients
        W=ones(1,N);

        %distance between all nodes (shortest paths)
        D_full=distances(graph(full_adjacency_mat));
        D_full(D_full==Inf)=10^10;
        D=D_full(1:N,1:N);
        D2=distances(graph(adjacency_mat));
        D2(D2==Inf)=10^10;
        
        %compute reduction in betas due to transcendence across all nodes
        try
            transcending_immunity= 1- exp(-D(vaccine_vector,:)/transcendence)';
        catch
            vaccine_vector;
            disp('ERROR')
        end
        %weight effects by vaccine strength (from 0 to 1 for each vaccine
        % strain).
        %transcending_immunity= 1- vaccine_vector(find(vaccine_vector)) .* (1-transcending_immunity);

        %take product for overlapping effects of multiple vaccination strains
        W=prod(transcending_immunity,2)';

        % THE FITNESS:
        %find proportion of supercritical strains: minimize this proportion
        node_remove=find(W<threshold);
        G=graph(adjacency_mat);
        G=rmnode(G,node_remove);
        [~, binsize]=conncomp(G);
        mean_compsize_P(loop)=sum(binsize.^2)/(N^2);
        
    end
    
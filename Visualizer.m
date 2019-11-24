function mean_compsize_P=Visualizer(vaccine_vector, adjacency_mat, threshold, transcendence)
    
    
    % # nodes
    N=size(adjacency_mat,1);

    %strain transmission coefficients
    W=ones(1,N);

    %distance between all nodes (shortest paths)
    D=distances(graph(adjacency_mat));

    %compute reduction in betas due to transcendence across all nodes
    try
        transcending_immunity= 1- exp(-D(vaccine_vector,:)/transcendence)';
    catch
        vaccine_vector;
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
        
    %PLOT
    figure; 
    H1=graph(adjacency_mat);
    %H2=rmnode(H1,node_remove);
    h=plot(H1,'layout','force','MarkerSize',10,'NodeLabel',{},'NodeColor','k',...
        'EdgeColor','k','LineWidth',5);

    H2=H1;
    for i=1:length(node_remove)
        node_idx=node_remove(i);
        nbors = neighbors(H1,node_idx);
        highlight(h,node_idx,nbors,'EdgeColor','r','LineWidth',1)
        H2=H2.rmedge(node_idx,nbors);
    end
    highlight(h,node_remove,'NodeColor','r','MarkerSize',5);
    highlight(h,vaccine_vector,'NodeColor','b','MarkerSize',15);
    
    
% % %     subplot(1,2,2)
% % %     G=rmnode(G,node_remove);
% % %     [~, binsize]=conncomp(G);
% % %     g=plot(H2,'layout','force','MarkerSize',10,'UseGravity',true,'NodeLabel',{});
% % %     highlight(g,node_remove,'NodeColor','r','MarkerSize',5);
% % %     highlight(g,vaccine_vector,'NodeColor','b','MarkerSize',15);
% % %     %highlight(g,node_remove,'NodeColor','k','MarkerSize',5);
% % %     %highlight(g,vaccine_vector,'NodeColor','r','MarkerSize',15);
% % %     p_fit=sum(binsize.^2)/(N^2);
% % %     title(p_fit);
 
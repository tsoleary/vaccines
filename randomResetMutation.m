function mutantChildren = randomResetMutation(parents,options,GenomeLength,FitnessFcn,state,thisScore,population,probMutating)

    mutantChildren = population(parents,:);
    
    for i=1:size(parents)
        for j=1:size(population,2)
       
            %mutate
            if rand()<probMutating
                N=1/probMutating;
                mutantChildren(i,j)=ceil(rand()*(N));
% % %                 if rand()<1
% % %                    mutantChildren(i,j)=ceil(rand()*(N));
% % %                 else
% % %                     try
% % %                         cur=mutantChildren(i,j);
% % %                         nbors=find(current_adjacency(cur,:));
% % %                         mutantChildren(i,j)=datasample(nbors,1);
% % %                     catch
% % %                         continue
% % %                     end
% % %                 end
            end
        end
    end
end
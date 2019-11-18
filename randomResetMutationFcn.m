function mutantChildren = randomResetMutationFcn(parents,options,GenomeLength,FitnessFcn,state,thisScore,population,probMutating)

    mutantChildren = population(parents,:);
    
    for i=1:size(parents)
        for j=1:size(population,2)
       
            %mutate
            if rand()<probMutating
               mutantChildren(i,j)=ceil(rand()*(1/probMutating));
            end
        end
    end
end
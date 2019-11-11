function normGenom = NormalizeGenome(genomeMat, thresh) 

% NormalizeGenome takes a binary genome and a threashold value and returns
% a normalized genome that sums to the threshold value. If the number of
% ones in a genome is less than the threshold value, the genome is left
% untouched to avoid alleles that excede 1

% INTPUTS: 
%       genomeMat: a binary matrix rows = genomes, cols = individuals
%       thresh: a scaler value that each genome should sum to
% OUTPUTS: 
%       normGenom: a real number matrix (0 to 1) rows = genomes, cols = individuals
%
% Written onNovember 10th 2019 by Alex Burnham

    popSize = length(genomeMat(:,1)); % calc the pop size
    genSize = length(genomeMat(1,:)); % calc the genome size

    normGenom = zeros(popSize, genSize); % init normalized genome

for i = 1:popSize
    
    % if the genom has greater than 3 or more ones, normalize
    if sum(genomeMat(i,:)) >= thresh
    normGenom(i,:) = (genomeMat(i,:)./sum(genomeMat(i,:))).*thresh;
    
    % if it doesn't, leave it alone to avoid being outsize a 0 to 1 boundery 
    else
    normGenom(i,:) = (genomeMat(i,:)); 
    
    end
end

end

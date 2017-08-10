
function tempPerClu = findTempForEachClu(clu, spikeTemplates)
% function tempPerClu = findTempForEachClu(clu, spikeTemplates)
%
% determine which template "corresponds to" each cluster, meaning, which
% template is most represented for each cluster
%
% output:
% - tempPerClu: a vector such that tempPerClu(clusterID) = templateID,
% where templateID is the template most represented for that cluster. E.g.
% if tempPerClu(1001) = 15, it means that the template most represented
% among the spikes for clu==1000 is template 15. templateID will be
% zero-indexed (as they were in spikeTemplates), but the entries of 
% tempPerClu cannot be (so they are shifted by one), which is why cluster 
% 1000 is found at index 1001. 

tempCountsByClu = full(sparse(double(clu)+1, double(spikeTemplates)+1, ones(size(clu))));

[~, tempPerClu] = max(tempCountsByClu,[],2);
tempPerClu = tempPerClu-1;
tempPerClu(~ismember(0:numel(tempPerClu)-1, unique(clu))) = NaN; % these entries will all have been "0" but they are clusters that don't exist

% old algorithm, was not realizing "max" will take care of this
% tempCountsByClu = tempCountsByClu+rand(size(tempCountsByClu))/10; % now there's just one max
% [tempPerClu, ~] = find(bsxfun(@eq, tempCountsByClu, max(tempCountsByClu,[],2))');
% tempPerClu = tempPerClu-1;

function tempPerClu = findTempForEachClu(clu, spikeTemplates)
% function tempPerClu = findTempForEachClu(clu, spikeTemplates)
%
% determine which template "corresponds to" each cluster, meaning, which
% template is most represented for each cluster

tempCountsByClu = full(sparse(double(clu)+1, double(spikeTemplates)+1, ones(size(clu))));
tempCountsByClu = tempCountsByClu+rand(size(tempCountsByClu))/10; % now there's just one max
[tempPerClu, ~] = find(bsxfun(@eq, tempCountsByClu, max(tempCountsByClu,[],2))');
tempPerClu = tempPerClu-1;
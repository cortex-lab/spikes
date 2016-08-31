

function [spikeWidths, tempWidths, clusterWidths] = computeSpikeWidths(tempsUnW, spikeTemplates, varargin)
% function [spikeWidths, tempWidths, clusterWidths] = computeSpikeWidths(tempsUnW, spikeTemplates[, clu])
% Computes spike widths in samples, according to widths of each template
% If you pass clu, will also compute it for clusters.


% The amplitude on each channel is the positive peak minus the negative
tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));

% The template amplitude is the amplitude of its largest channel
% tempAmps = max(tempChanAmps,[],2);

nTemps = size(tempChanAmps,1);
tempWidths = zeros(1, nTemps);
for t = 1:nTemps
    theseChanAmps = tempChanAmps(t,:);
    maxWF = tempsUnW(t,:,find(theseChanAmps==max(theseChanAmps),1));
    tempWidths(t) = find(maxWF==max(maxWF),1)-find(maxWF==min(maxWF),1);
end

spikeWidths = tempWidths(spikeTemplates+1)';

if ~isempty(varargin)
    clu = varargin{1};
    clusterWidths = clusterAverage(clu, spikeWidths);
else
    clusterWidths = [];
end
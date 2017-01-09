

function [timeBins, depthBins, allP] = psthByDepth(spikeTimes, spikeDepths, depthBinSize, timeBinSize, eventTimes, win, bslWin)
% function [timeBins, depthBins, allP] = psthByDepth(spikeTimes, spikeDepths, depthBinSize, timeBinSize, eventTimes, win, bslWin)
%
% Computes PSTH split by the depths of spikes provided
% If pass a "bslWin" will normalize to the bin counts in that period
% - if so normalized, units of allP is stdev relative to baseline mean
% - if not normalized, units of allP is spikes/sec
%
% timeBins is 1xnTimeBins
% depthBins is 1xnDepthBins
% allP is nDepthBins x nTimeBins
% 


depthBins = min(spikeDepths):depthBinSize:max(spikeDepths); nD = length(depthBins)-1;
bslEventTimes = eventTimes; % could make it so you can provide a different bsl event

for d = 1:nD
    theseSp = spikeDepths>depthBins(d) & spikeDepths<=depthBins(d+1);
    
    if ~isempty(bslWin)
        [psth, ~, ~, ~, ~, ~] = psthAndBA(spikeTimes(theseSp), bslEventTimes, bslWin, timeBinSize);
        normMn = mean(psth); normStd = std(psth);
    end
    
    [psth, timeBins, ~, ~, ~, ~] = psthAndBA(spikeTimes(theseSp), eventTimes, win, timeBinSize);
    
    if d==1
        allP = zeros(nD, length(psth));
    end
    if ~isempty(bslWin) && normStd>0
        allP(d,:) = (psth-normMn)./normStd;
    else
        allP(d,:) = psth;
    end
        
end

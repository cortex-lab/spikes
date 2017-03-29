

function [timeBins, depthBins, allP, normVals] = psthByDepth(spikeTimes, spikeDepths, depthBinSize, timeBinSize, eventTimes, win, bslWin, varargin)
% function [timeBins, depthBins, allP] = psthByDepth(spikeTimes, ...
%   spikeDepths, depthBinSize, timeBinSize, eventTimes, win, bslWin[, bslEvents])
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
if ~isempty(varargin)
    bslEventTimes = varargin{1};
else
    bslEventTimes = eventTimes;
end

if ~isempty(bslWin)
    normVals = zeros(nD,2);
else
    normVals = [];
end

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
        normVals(d,:) = [normMn normStd];
    else
        allP(d,:) = psth;
    end
        
end

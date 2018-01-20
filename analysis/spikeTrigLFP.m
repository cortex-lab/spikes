

function mnLFP = spikeTrigLFP(tLFP, lfpdat, theseST, winAroundSpike)
% function mnLFP = spikeTrigLFP(tLFP, lfpdat, theseST, winAroundSpike)
%
% returns nChannels x nTimePoints mean spike-triggered LFP. 
%
% Inputs:
% - tLFP - 1 x nLFP - vector of time points at which the lfp was sampled
% (s)
% - lfpdat - nChannels x nLFP - the lfp data (any units)
% - theseST - nSpikes x 1 - spike times (s)
% - winAroundSpike - 1 x nTimePoints - time points around each spike to
% include

% the time points to sample LFP at
sampTimes = bsxfun(@plus, theseST, winAroundSpike);
    
mnLFP = zeros(size(lfpdat,1), numel(winAroundSpike));

% compute the spike-trig LFP channel-by-channel
%     for ch = 1:size(lfpdat,1)
%         if mod(ch,10)==0
%             fprintf(1, 'ch %d/%d...\n', ch, size(lfpdat,1));
%         end
%         mnLFP(ch,:) = mean(interp1(tLFP, lfpdat(ch,:)', sampTimes),1);
%     end

% compute the spike-trig LFP chunk-by-chunk (this one seems faster)
nChunk = 4;
chunkBounds = linspace(0,max(tLFP),nChunk+1);
for n = 1:nChunk
    fprintf(1, 'chunk %d/%d...\n', n, nChunk);
    inclSamps = tLFP>chunkBounds(n) & tLFP<=chunkBounds(n+1);
    inclSpikes = theseST>chunkBounds(n) & theseST<=chunkBounds(n+1);
    mnLFP = mnLFP+squeeze(nanmean(interp1(tLFP(inclSamps), ...
        lfpdat(:,inclSamps)', sampTimes(inclSpikes,:)),1))';
end
mnLFP = mnLFP/nChunk;
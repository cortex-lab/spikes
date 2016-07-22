
function [pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths, ampBins, depthBins, recordingDur);

nDbins = length(depthBins)-1;
nAbins = length(ampBins)-1;

pdfs = zeros(nDbins, nAbins);
cdfs = zeros(nDbins, nAbins);
for b = 1:length(depthBins)-1
    h = histc(spikeAmps(spikeDepths>depthBins(b) & spikeDepths<=depthBins(b+1)), ampBins);
    h = h./recordingDur; % convert to Hz
    pdfs(b,:) = h(1:end-1);
    thiscdf = cumsum(h(end-1:-1:1));
    cdfs(b,:) = thiscdf(end:-1:1);
end
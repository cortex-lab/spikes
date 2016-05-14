


function [binArray, binCenters] = timestampsToBinned(timeStamps, referencePoints, binSize, window)
% [binArray, binCenters] = timestampsToBinned(timeStamps, referencePoints, binSize,
% window)
%
% Returns an array of binned spike counts. If you use a large enough
% binSize, it may well be possible that there is more than one spike in a
% bin, i.e. the value of some bin may be >1. 

binBorders = window(1):binSize:window(2);
numBins = length(binBorders)-1;

if isempty(referencePoints)
    binArray = [];
    binCenters = binBorders(1:end-1)+binSize/2; % not sure if this is the same as what you get below?
    return;
end

binArray = zeros(length(referencePoints), numBins);

if isempty(timeStamps)
    binCenters = binBorders(1:end-1)+binSize/2; % not sure if this is the same as what you get below?
    return;
end

for r = 1:length(referencePoints)
    [n,binCenters] = histdiff(timeStamps, referencePoints(r), binBorders);
    binArray(r,:) = n;
end

    


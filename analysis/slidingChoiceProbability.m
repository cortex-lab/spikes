

function [cp, t, rs] = slidingChoiceProbability(ba, baBinSize, cpBinSize, trLabels)
% function [cp, t, rs] = slidingChoiceProbability(ba, baBinSize, cpBinSize, trLabels)
% Choice probability (i.e. area under ROC curve between spike counts in a
% bin) using a sliding window. Input argument "ba" is a "binnedArray" from
% timestampsToBinned or psthAndBA. baBinSize is the bin size of the columns
% of ba; cpBinSize is the desired bin size for counting spikes. trLabels
% must only have two different elements. They will be ordered as returned
% by "unique(trLabels)", and if the first has larger mean, then cp<0.5. If
% a third output argument is requested, it is the ranksum p-value.

% time indices of the ba
baT = (0:size(ba,2)-1)*baBinSize;

uLabels = unique(trLabels);
assert(length(uLabels)==2, 'For slidingChoiceProbability, exactly two labels required');

tr1 = trLabels==uLabels(1);
tr2 = trLabels==uLabels(2);

% compute counts with box convolution
nBins = round(cpBinSize/baBinSize);
bx = ones(1, nBins)/nBins; % box filter with area=1 so units don't change
counts = conv2(1, bx, ba, 'same');
counts1 = counts(tr1,:);
counts2 = counts(tr2,:);

% time indices to compute the cp for
tIncl = find(baT>cpBinSize/2 & baT<max(baT)-cpBinSize/2);
t = baT(tIncl);
cp = zeros(1, length(t));
rs = zeros(1, length(t));
for tInd = 1:length(t)

    c1 = counts1(:,tIncl(tInd));
    c2 = counts2(:,tIncl(tInd));

    cp(tInd) = rocArea(c1, c2);
    
    if nargout>2
        rs(tInd) = ranksum(c1, c2);
    end
end

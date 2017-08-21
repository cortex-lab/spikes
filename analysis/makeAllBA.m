

function [allba, bins] = makeAllBA(st, clu, eventTimes, binSize, win)

cIDs = unique(clu);

binBorders = win(1):binSize:win(2);
numBins = length(binBorders)-1;

allba = zeros(numel(cIDs), numel(eventTimes), numBins);
% there will be a faster way to implement this with WithinRanges, at least
% in the case of non-overlapping ranges
for c = 1:numel(cIDs)
%     fprintf(1, '%d/%d\n', c, numel(cIDs))
    theseST = st(clu==cIDs(c));

    [ba, bins] = timestampsToBinned(theseST, eventTimes, binSize, win);
    
    allba(c, :,:) = ba;
end